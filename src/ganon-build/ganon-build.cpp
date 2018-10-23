#include <utils/safequeue.hpp>

#include <cxxopts.hpp>
#include <seqan/kmer.h>

#include <chrono>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <vector>
#include <tuple>

static const uint64_t gbInBits = 8589934592;

struct Seqs
{
    std::string acc;
    seqan::Dna5String seq;
};

int main( int argc, char* argv[] )
{
    int old_argc = argc; // parser always set argc to 1

    cxxopts::Options options( "ganon-build", "Ganon bloom filter build" );

    // clang-format off
    options.add_options()
        ( "e,seqid-bin", "Seqid bin file", cxxopts::value< std::string >() )
        ( "o,output-file", "Output file", cxxopts::value< std::string >() )
        ( "s,bloom-size", "Final bloom filter size in GB", cxxopts::value< int >()->default_value( "16" ) )
        ( "bloom-size-bits", "Final bloom filter size in bits", cxxopts::value< uint64_t >()->default_value( "0" ) )
        ( "k,kmer-size", "K size", cxxopts::value< int >()->default_value( "19" ) )
        ( "n,hash-functions", "Number of hash functions", cxxopts::value< int >()->default_value( "3" ) )
        ( "t,threads", "Number of threads", cxxopts::value< int >()->default_value( "1" ) )
        ( "h,help", "Print help" )( "v,version", "Show version" )
        ( "references", "references", cxxopts::value< std::vector< std::string > >() );
    // clang-format on

    options.parse_positional( { "references" } );
    options.positional_help( "reference.fna[.gz] [reference2.fna[.gz] ... referenceN.fna[.gz]]" );

    auto args = options.parse( argc, argv );

    if ( args.count( "help" ) || old_argc <= 1 )
    {
        std::cerr << options.help() << std::endl;
        return 0;
    }
    else if ( args.count( "version" ) )
    {
        std::cerr << "version" << std::endl;
        return 0;
    }

    uint64_t bloom_filter_size;
    if ( args["bloom-size-bits"].as< uint64_t >() > 0 )
    {
        bloom_filter_size = args["bloom-size-bits"].as< uint64_t >();
    }
    else
    {
        bloom_filter_size = args["bloom-size"].as< int >() * gbInBits; // gbInBits -> 8589934592 bits = 1 Gb
    }

    std::cerr << "seqid-bin: " << args["seqid-bin"].as< std::string >() << std::endl;
    std::cerr << "bloom-size: " << std::fixed << std::setprecision( 2 ) << (float) bloom_filter_size / (float) gbInBits
              << std::endl;
    std::cerr << "bloom-size-bits: " << bloom_filter_size << std::endl;
    std::cerr << "kmer-size: " << args["kmer-size"].as< int >() << std::endl;
    std::cerr << "hash-functions: " << args["hash-functions"].as< int >() << std::endl;
    std::cerr << "threads: " << args["threads"].as< int >() << std::endl;
    std::cerr << "output-file: " << args["output-file"].as< std::string >() << std::endl;
    std::cerr << "references: " << std::endl;
    for ( const auto& s : args["references"].as< std::vector< std::string > >() )
    {
        std::cerr << s << std::endl;
    }

    uint64_t kmer_size      = args["kmer-size"].as< int >();
    int      hash_functions = args["hash-functions"].as< int >();

    int threads = args["threads"].as< int >();

    std::map< std::string, std::vector< std::tuple<uint64_t,uint64_t,uint64_t> > > bins;
    std::ifstream                     infile( args["seqid-bin"].as< std::string >() );
    std::string                       seqid;
    uint64_t                          seqstart;
    uint64_t                          seqend;
    uint64_t                          bin;
    uint64_t                          noBins = 0;
    while ( infile >> seqid >> seqstart >> seqend >> bin )
    {
        bins[seqid].push_back(std::make_tuple(seqstart, seqend, bin));
        if ( bin > noBins )
            noBins = bin;
    }
    noBins = noBins + 1;
    std::cerr << bins.size() << " sequences on " << noBins << " bins" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    seqan::KmerFilter< seqan::Dna5, seqan::InterleavedBloomFilter, seqan::Uncompressed > filter(
        noBins, hash_functions, kmer_size, bloom_filter_size );

    std::chrono::duration< double > elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Creating Bloom filter: " << elapsed.count() << std::endl;

    std::mutex                         mtx;
    SafeQueue< Seqs >                  q;
    std::vector< std::future< void > > tasks;
    bool                               finished = false;

    // Start threads async
    start = std::chrono::high_resolution_clock::now();
    for ( int taskNo = 0; taskNo < threads; ++taskNo )
    {
        tasks.emplace_back( std::async( std::launch::async, [=, &bins, &filter, &q, &finished, &mtx] {
            while ( true )
            {
                Seqs val = q.pop();
                if ( val.acc != "" )
                { // if not empty
                    for (uint64_t i=0; i < bins[val.acc].size(); i++)
                    {
                        auto [ fragstart, fragend, binid ] = bins[val.acc][i];
                        // For infixes, we have to provide both the including start and the excluding end position.
                        // fragstart -1 to fix offset
                        // fragend -1+1 to fix offset and not exclude last position
                        seqan::Infix< seqan::Dna5String >::Type fragment = infix(val.seq, fragstart - 1, fragend);
                        seqan::insertKmer( filter, fragment, binid, 0 );
                        mtx.lock();
                        std::cerr << val.acc << " [" << fragstart << ":" << fragend << "] added to bin " << binid << std::endl;
                        std::cerr << fragment << std::endl;
                        mtx.unlock();
                    }
                }
                if ( finished && q.empty() )
                    break;
            }
        } ) );
    }

    // extra thread for reading the input
    tasks.emplace_back( std::async( std::launch::async, [=, &bins, &q, &finished, &mtx] {
        for ( auto const& reference_fasta_file : args["references"].as< std::vector< std::string > >() )
        {
            seqan::SeqFileIn seqFileIn;
            if ( !seqan::open( seqFileIn, seqan::toCString( reference_fasta_file ) ) )
            {
                std::cerr << "Unable to open " << reference_fasta_file << std::endl;
                continue;
            }
            while ( !seqan::atEnd( seqFileIn ) )
            {
                while ( q.size() > ( threads * 10 ) )
                {
                    ; // spin
                }
                seqan::StringSet< seqan::CharString >  ids;
                seqan::StringSet< seqan::IupacString > seqs;
                seqan::readRecords( ids, seqs, seqFileIn, threads * 5 );
                for ( uint16_t seqID = 0; seqID < seqan::length( ids ); ++seqID )
                {
                    if ( seqan::length( seqs[seqID] ) < kmer_size )
                    { // sequence too small
                        mtx.lock();
                        std::cerr << ids[seqID] << " has sequence smaller than k-mer size" << std::endl;
                        mtx.unlock();
                        continue;
                    }
                    std::string cid = seqan::toCString( ids[seqID] );
                    std::string acc = cid.substr( 0, cid.find( ' ' ) );
                    if ( bins.count( acc ) == 0 )
                    { // not defined on the bins
                        mtx.lock();
                        std::cerr << acc << " not defined on bins file" << std::endl;
                        mtx.unlock();
                        continue;
                    }
                    q.push( Seqs{ acc, seqs[seqID] } );
                }
            }
            seqan::close( seqFileIn );
        }
        finished = true;
    } ) );

    for ( auto&& task : tasks )
    {
        task.get();
    }
    elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Adding k-mers: " << elapsed.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    seqan::store( filter, seqan::toCString( args["output-file"].as< std::string >() ) );
    elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Saving Bloom filter: " << elapsed.count() << std::endl;

    return 0;
}
