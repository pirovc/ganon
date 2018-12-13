#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <seqan/kmer.h>
#include <tuple>
#include <utils/safequeue.hpp>
#include <vector>

#include "Arguments.hpp"

struct Seqs
{
    std::string       seqid;
    seqan::Dna5String seq;
};


struct FragmentBin
{
    uint64_t start;
    uint64_t end;
    uint64_t bin_id;
};

typedef std::map< std::string, std::vector< FragmentBin > >                                  TSeqBin;
typedef seqan::KmerFilter< seqan::Dna5, seqan::InterleavedBloomFilter, seqan::Uncompressed > TInterleavedBloomFilter;

void parse_seqid_bin( const std::string& seqid_bin_file, TSeqBin& seq_bin, std::set< uint64_t >& bin_ids )
{
    std::string   line;
    std::ifstream infile( seqid_bin_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field );
        // seqid <tab> seqstart <tab> seqend <tab> binid
        uint32_t binid = std::stoul( fields[3] );
        seq_bin[fields[0]].push_back( FragmentBin{ std::stoul( fields[1] ), std::stoul( fields[2] ), binid } );
        bin_ids.insert( binid );
    }
}


TInterleavedBloomFilter load_filter( Arguments& args, const std::set< uint64_t >& bin_ids )
{

    uint64_t number_of_bins;
    if ( !args.update_filter_file.empty() )
    {
        // dummy variables in case of update
        number_of_bins = args.hash_functions = args.kmer_size = args.filter_size = 1;
    }
    else
    {
        number_of_bins = bin_ids.size();
    }

    // define filter
    TInterleavedBloomFilter filter( number_of_bins, args.hash_functions, args.kmer_size, args.filter_size );

    // load from disk in case of update
    if ( !args.update_filter_file.empty() )
    {
        // load filter
        seqan::retrieve( filter, seqan::toCString( args.update_filter_file ) );

        args.kmer_size = seqan::getKmerSize( filter );

        uint32_t number_of_bins_before = seqan::getNumberOfBins( filter );

        // Reset bins if complete set of sequences is provided (re-create updated bins)
        if ( args.update_complete )
        {
            std::vector< uint32_t > updated_bins;
            updated_bins.insert( updated_bins.end(), bin_ids.begin(), bin_ids.end() );
            seqan::clear( filter, updated_bins, args.threads );
        }

        // TODO -> create new bins on the loaded filter
        std::cerr << *bin_ids.rbegin() + 1 - number_of_bins_before << " new bins added to existing  "
                  << number_of_bins_before << " bins" << std::endl;
    }
    return filter;
}

int main( int argc, char** argv )
{

    // Parse arguments
    Arguments args( argc, argv );
    if ( !args.parse() )
        return 0;

    if ( args.verbose )
        args.print();
    //////////////////////////////


    // Load seqid-bin and filter
    auto start = std::chrono::high_resolution_clock::now();

    // parse seqid bin
    TSeqBin              seq_bin;
    std::set< uint64_t > bin_ids;
    parse_seqid_bin( args.seqid_bin_file, seq_bin, bin_ids );
    std::cerr << seq_bin.size() << " sequences, " << bin_ids.size() << " bins defined on " << args.seqid_bin_file
              << std::endl;

    // load new or given filter
    TInterleavedBloomFilter filter = load_filter( args, bin_ids );

    std::chrono::duration< double > elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Loading filter and files: " << elapsed.count() << "s" << std::endl;
    //////////////////////////////


    // Start execution threads to add kmers
    std::mutex                         mtx;
    SafeQueue< Seqs >                  q_Seqs;
    std::vector< std::future< void > > tasks;
    bool                               finished = false;
    start                                       = std::chrono::high_resolution_clock::now();
    for ( int taskNo = 0; taskNo < args.threads; ++taskNo )
    {
        tasks.emplace_back( std::async( std::launch::async, [=, &seq_bin, &filter, &q_Seqs, &finished, &mtx, &args] {
            while ( true )
            {
                Seqs val = q_Seqs.pop();
                if ( val.seqid != "" )
                {
                    for ( uint64_t i = 0; i < seq_bin[val.seqid].size(); i++ )
                    {
                        auto [fragstart, fragend, binid] = seq_bin[val.seqid][i];
                        // For infixes, we have to provide both the including start and the excluding end position.
                        // fragstart -1 to fix offset
                        // fragend -1+1 to fix offset and not exclude last position
                        seqan::Infix< seqan::Dna5String >::Type fragment = infix( val.seq, fragstart - 1, fragend );
                        seqan::insertKmer( filter, fragment, binid, 0 );
                        if ( args.verbose )
                        {
                            mtx.lock();
                            std::cerr << val.seqid << " [" << fragstart << ":" << fragend << "] added to bin " << binid
                                      << std::endl;
                            mtx.unlock();
                        }
                    }
                }
                if ( finished && q_Seqs.empty() )
                    break;
            }
        } ) );
    }
    //////////////////////////////

    // Start extra thread for reading the input
    tasks.emplace_back( std::async( std::launch::async, [=, &seq_bin, &q_Seqs, &finished, &mtx] {
        for ( auto const& reference_file : args.reference_files )
        {
            seqan::SeqFileIn seqFileIn;
            if ( !seqan::open( seqFileIn, seqan::toCString( reference_file ) ) )
            {
                std::cerr << "Unable to open " << reference_file << std::endl;
                continue;
            }
            while ( !seqan::atEnd( seqFileIn ) )
            {
                while ( q_Seqs.size() > ( args.threads * 10 ) )
                {
                    ; // spin
                }
                seqan::StringSet< seqan::CharString >  ids;
                seqan::StringSet< seqan::IupacString > seqs;
                seqan::readRecords( ids, seqs, seqFileIn, args.threads * 5 );
                for ( uint64_t i = 0; i < seqan::length( ids ); ++i )
                {
                    if ( seqan::length( seqs[i] ) < args.kmer_size )
                    { // sequence too small
                        mtx.lock();
                        std::cerr << ids[i] << " has sequence smaller than k-mer size" << std::endl;
                        mtx.unlock();
                        continue;
                    }
                    std::string cid   = seqan::toCString( ids[i] );
                    std::string seqid = cid.substr( 0, cid.find( ' ' ) );
                    if ( seq_bin.count( seqid ) == 0 )
                    {
                        mtx.lock();
                        std::cerr << seqid << " not defined on seqid-bin file" << std::endl;
                        mtx.unlock();
                        continue;
                    }
                    q_Seqs.push( Seqs{ seqid, seqs[i] } );
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
    std::cerr << "Adding k-mers: " << elapsed.count() << "s" << std::endl;
    //////////////////////////////

    // Store filter
    start = std::chrono::high_resolution_clock::now();
    seqan::store( filter, seqan::toCString( args.output_filter_file ) );
    elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Saving filter: " << elapsed.count() << "s" << std::endl;
    //////////////////////////////

    return 0;
}
