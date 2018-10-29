#include <defaults/defaults.hpp>
#include <utils/SafeQueue.hpp>

#include <cxxopts.hpp>
#include <seqan/kmer.h>

#include <chrono>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <tuple>
#include <vector>

static const uint64_t gbInBits = 8589934592;

struct Seqs
{
    std::string       acc;
    seqan::Dna5String seq;
};

int main( int argc, char* argv[] )
{
    int old_argc = argc; // parser always set argc to 1

    cxxopts::Options options( "ganon-build", "Ganon bloom filter build" );

    // clang-format off
    options.add_options()
        ( "e,seqid-bin", "Space-separated file with the following fields: Seq. Identifier <space> Pos. Start <space> Pos. End <space> Bin Id", cxxopts::value< std::string >() )
        ( "o,output-file", "Output file", cxxopts::value< std::string >() )
        ( "s,bloom-size", "Final bloom filter size in GB", cxxopts::value< uint64_t >()->default_value( "16" ) )
        ( "b,bloom-size-bits", "Final bloom filter size in bits", cxxopts::value< uint64_t >()->default_value( "0" ) )
        ( "k,kmer-size", "k size", cxxopts::value< uint64_t >()->default_value( "19" ) )
        ( "n,hash-functions", "Number of hash functions", cxxopts::value< uint16_t >()->default_value( "3" ) )
        ( "u,update-bloom-filter", "If provided, filte updated with new sequences", cxxopts::value< std::string >()->default_value( "" ) )
        ( "c,update-complete", "Old and new sequences are provided for updated bins", cxxopts::value< bool >()->default_value( "false" ) )
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >()->default_value( "1" ) )
        //( "silent", "Silent mode", cxxopts::value<bool>()->default_value("false"))
        //( "verbose", "Verbose mode to STDERR", cxxopts::value<bool>()->default_value("false"))
        ( "h,help", "Print help" )
        ( "v,version", "Show version" )
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
        std::cerr << "version: " << defaults::version_string << std::endl;
        return 0;
    }

    std::string seqid_bin_file           = args["seqid-bin"].as< std::string >();
    uint64_t    bloom_filter_size        = 0;
    uint64_t    kmer_size                = args["kmer-size"].as< uint64_t >();
    uint16_t    hash_functions           = args["hash-functions"].as< uint16_t >();
    std::string update_bloom_filter_file = args["update-bloom-filter"].as< std::string >();
    bool        update_complete          = args["update-complete"].as< bool >();
    uint16_t    threads                  = args["threads"].as< uint16_t >();
    std::string output_file              = args["output-file"].as< std::string >();


    uint64_t number_of_bins;

    // Updating
    if ( !update_bloom_filter_file.empty() )
    {
        std::cerr
            << "--bloom-size[-bits], --kmer-size --hash-funtions ignored, using metadata from --update-bloom-filter"
            << std::endl;
    }
    else
    {
        if ( args["bloom-size-bits"].as< uint64_t >() > 0 )
        {
            bloom_filter_size = args["bloom-size-bits"].as< uint64_t >();
        }
        else
        {
            bloom_filter_size = args["bloom-size"].as< uint64_t >() * gbInBits; // gbInBits -> 8589934592 bits = 1 Gb
        }
    }

    std::cerr << "seqid-bin: " << seqid_bin_file << std::endl;
    std::cerr << "bloom-size: " << std::fixed << std::setprecision( 2 ) << (float) bloom_filter_size / (float) gbInBits
              << std::endl;
    std::cerr << "bloom-size-bits: " << bloom_filter_size << std::endl;
    std::cerr << "kmer-size: " << kmer_size << std::endl;
    std::cerr << "hash-functions: " << hash_functions << std::endl;
    std::cerr << "bloom-filter: " << update_bloom_filter_file << std::endl;
    std::cerr << "update-complete: " << update_complete << std::endl;
    std::cerr << "threads: " << threads << std::endl;
    std::cerr << "output-file: " << output_file << std::endl;
    std::cerr << "references: " << std::endl;
    for ( const auto& s : args["references"].as< std::vector< std::string > >() )
    {
        std::cerr << s << std::endl;
    }

    std::map< std::string, std::vector< std::tuple< uint64_t, uint64_t, uint64_t > > > bins;
    std::ifstream                                                                      infile( seqid_bin_file );
    std::string                                                                        seqid;
    uint64_t                                                                           seqstart;
    uint64_t                                                                           seqend;
    uint64_t                                                                           bin;
    std::unordered_set< uint64_t >                                                     updated_bins;
    uint32_t                                                                           max_bin_updated = 0;

    while ( infile >> seqid >> seqstart >> seqend >> bin )
    {
        bins[seqid].push_back( std::make_tuple( seqstart, seqend, bin ) );
        updated_bins.insert( bin );
        if ( bin > max_bin_updated )
        {
            max_bin_updated = bin;
        }
    }
    std::cerr << bins.size() << " sequences will be added on " << updated_bins.size() << " bins" << std::endl;
    number_of_bins = updated_bins.size();

    auto start = std::chrono::high_resolution_clock::now();
    typedef seqan::KmerFilter< seqan::Dna5, seqan::InterleavedBloomFilter, seqan::Uncompressed > IBloomFilter;

    // dummy variables in case of update
    if ( !update_bloom_filter_file.empty() )
    {
        number_of_bins = hash_functions = kmer_size = bloom_filter_size = 1;
    }


    // define filter
    IBloomFilter filter( number_of_bins, hash_functions, kmer_size, bloom_filter_size );

    // load from disk in case of update
    if ( !update_bloom_filter_file.empty() )
    {
        // load filter
        retrieve( filter, seqan::toCString( update_bloom_filter_file ) );
        kmer_size = seqan::getKmerSize( filter );

        uint32_t number_of_bins_before = seqan::getNumberOfBins( filter );

        // Reset bins if complete set of sequences is provided (re-create updated bins)
        if ( update_complete )
        {
            std::vector< uint32_t > ubins;
            ubins.insert( ubins.end(), updated_bins.begin(), updated_bins.end() );
            seqan::clear( filter, ubins, threads );
        }

        // TODO -> create new bins on the loaded bloom
        std::cerr << max_bin_updated + 1 - number_of_bins_before << " new bins added to existing  "
                  << number_of_bins_before << " bins" << std::endl;
    }

    std::chrono::duration< double > elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Creating/Loading Bloom filter: " << elapsed.count() << std::endl;

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
                    for ( uint64_t i = 0; i < bins[val.acc].size(); i++ )
                    {
                        auto [fragstart, fragend, binid] = bins[val.acc][i];
                        // For infixes, we have to provide both the including start and the excluding end position.
                        // fragstart -1 to fix offset
                        // fragend -1+1 to fix offset and not exclude last position
                        seqan::Infix< seqan::Dna5String >::Type fragment = infix( val.seq, fragstart - 1, fragend );
                        seqan::insertKmer( filter, fragment, binid, 0 );
                        mtx.lock();
                        std::cerr << val.acc << " [" << fragstart << ":" << fragend << "] added to bin " << binid
                                  << std::endl;
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
    seqan::store( filter, seqan::toCString( output_file ) );
    elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Saving Bloom filter: " << elapsed.count() << std::endl;

    return 0;
}
