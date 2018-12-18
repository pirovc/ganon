#include "GanonBuild.hpp"

#include <utils/SafeQueue.hpp>
#include <utils/Time.hpp>

#include <seqan/kmer.h>

#include <cinttypes>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace detail
{

struct Seqs
{
    std::string       seqid;
    seqan::Dna5String seq;
};

struct Stats
{
    Stats()
    : sumSeqLen{ 0 }
    , totalSeqsBinId{ 0 }
    , totalBinsBinId{ 0 }
    , totalSeqsFile{ 0 }
    , totalBinsFile{ 0 }
    , invalidSeqs{ 0 }
    , newBins{ 0 }
    {
    }

    uint64_t sumSeqLen;
    uint64_t totalSeqsBinId;
    uint32_t totalBinsBinId;
    uint64_t totalSeqsFile;
    uint32_t totalBinsFile;
    uint64_t invalidSeqs;
    uint32_t newBins;
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


TInterleavedBloomFilter load_filter( GanonBuild::Config& config, const std::set< uint64_t >& bin_ids, Stats& stats )
{
    uint64_t number_of_bins;
    if ( !config.update_filter_file.empty() )
    {
        // dummy variables in case of update
        number_of_bins = config.hash_functions = config.kmer_size = config.filter_size = 1;
    }
    else
    {
        number_of_bins = bin_ids.size();
    }

    // define filter
    TInterleavedBloomFilter filter( number_of_bins, config.hash_functions, config.kmer_size, config.filter_size );

    // load from disk in case of update
    if ( !config.update_filter_file.empty() )
    {
        // load filter
        seqan::retrieve( filter, seqan::toCString( config.update_filter_file ) );

        config.kmer_size = seqan::getKmerSize( filter );

        // Reset bins if complete set of sequences is provided (re-create updated bins)
        if ( config.update_complete )
        {
            std::vector< uint32_t > updated_bins;
            updated_bins.insert( updated_bins.end(), bin_ids.begin(), bin_ids.end() );
            seqan::clear( filter, updated_bins, config.threads );
        }

        // TODO -> create new bins on the loaded filter
        stats.newBins = *bin_ids.rbegin() + 1 - stats.totalBinsFile;
    }

    stats.totalBinsFile = seqan::getNumberOfBins( filter );

    return filter;
}

void print_time( const GanonBuild::Config& config,
                 Time&                     timeGanon,
                 Time&                     timeLoadFiles,
                 Time&                     timeLoadSeq,
                 Time&                     timeBuild,
                 Time&                     timeSaveFilter )
{
    std::cerr << "ganon-build       start time: " << timeGanon.get_start_ctime();
    std::cerr << "Loading files     start time: " << timeLoadFiles.get_start_ctime();
    std::cerr << "Loading files       end time: " << timeLoadFiles.get_end_ctime();
    std::cerr << "Loading sequences start time: " << timeLoadSeq.get_start_ctime();
    std::cerr << "Building          start time: " << timeBuild.get_start_ctime();
    std::cerr << "Loading sequences   end time: " << timeLoadSeq.get_end_ctime();
    std::cerr << "Building            end time: " << timeBuild.get_end_ctime();
    std::cerr << "Saving filter     start time: " << timeSaveFilter.get_start_ctime();
    std::cerr << "Saving filter       end time: " << timeSaveFilter.get_end_ctime();
    std::cerr << "ganon-build         end time: " << timeGanon.get_end_ctime();
    std::cerr << std::endl;
    std::cerr << " - loading files: " << timeLoadFiles.get_elapsed() << std::endl;
    std::cerr << " - loading sequences (1t): " << timeLoadSeq.get_elapsed() << std::endl;
    std::cerr << " - building filter (" << config.build_threads << "t): " << timeBuild.get_elapsed() << std::endl;
    std::cerr << " - saving filter: " << timeSaveFilter.get_elapsed() << std::endl;
    std::cerr << " - total: " << timeGanon.get_elapsed() << std::endl;
    std::cerr << std::endl;
}

void print_stats( Stats& stats, const GanonBuild::Config& config, Time& timeBuild )
{
    double   elapsed_build = timeBuild.get_elapsed();
    uint64_t validSeqs     = stats.totalSeqsFile - stats.invalidSeqs;
    std::cerr << "ganon-build processed " << validSeqs << " sequences (" << stats.sumSeqLen / 1000000.0 << " Mbp) in "
              << elapsed_build << " seconds (" << ( validSeqs / 1000.0 ) / ( elapsed_build / 60.0 ) << " Kseq/m, "
              << ( stats.sumSeqLen / 1000000.0 ) / ( elapsed_build / 60.0 ) << " Mbp/m)" << std::endl;
    std::cerr << " - " << stats.totalSeqsBinId << " sequences and " << stats.totalBinsBinId << " bins defined on "
              << config.seqid_bin_file << std::endl;
    std::cerr << " - " << stats.totalSeqsFile << " sequences (" << stats.invalidSeqs
              << " invalid) were read from the input sequence files." << std::endl;
    if ( !config.update_filter_file.empty() )
        std::cerr << " - " << stats.newBins << " new bins were added to the existing " << stats.totalBinsFile
                  << " bins." << std::endl;
    std::cerr << " - " << validSeqs << " valid sequences in " << stats.totalBinsFile + stats.newBins
              << " bins were written to " << config.output_filter_file << std::endl;
}

} // namespace detail

bool GanonBuild::run( Config config )
{
    Time timeGanon;
    timeGanon.start();
    Time timeLoadFiles;
    Time timeLoadSeq;
    Time timeBuild;
    Time timeSaveFilter;

    if ( config.verbose )
    {
        std::cerr << config;
    }

    //////////////////////////////

    detail::Stats stats;

    timeLoadFiles.start();
    // parse seqid bin
    detail::TSeqBin      seq_bin;
    std::set< uint64_t > bin_ids;
    parse_seqid_bin( config.seqid_bin_file, seq_bin, bin_ids );
    stats.totalSeqsBinId = seq_bin.size();
    stats.totalBinsBinId = bin_ids.size();

    // load new or given filter
    detail::TInterleavedBloomFilter filter = load_filter( config, bin_ids, stats );
    timeLoadFiles.end();
    //////////////////////////////

    // Start execution threads to add kmers
    timeBuild.start();
    std::mutex                         mtx;
    SafeQueue< detail::Seqs >          q_Seqs;
    std::vector< std::future< void > > tasks;
    bool                               finished = false;
    for ( int taskNo = 0; taskNo < config.build_threads; ++taskNo )
    {
        tasks.emplace_back( std::async( std::launch::async, [=, &seq_bin, &filter, &q_Seqs, &finished, &mtx, &config] {
            while ( true )
            {
                detail::Seqs val = q_Seqs.pop();
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
                        if ( config.verbose )
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
    timeLoadSeq.start();
    tasks.emplace_back( std::async( std::launch::async, [=, &seq_bin, &q_Seqs, &finished, &mtx, &timeLoadSeq, &stats] {
        for ( auto const& reference_file : config.reference_files )
        {
            seqan::SeqFileIn seqFileIn;
            if ( !seqan::open( seqFileIn, seqan::toCString( reference_file ) ) )
            {
                std::cerr << "Unable to open " << reference_file << std::endl;
                continue;
            }
            while ( !seqan::atEnd( seqFileIn ) )
            {
                while ( q_Seqs.size() > ( config.threads * 10 ) )
                {
                    ; // spin
                }
                seqan::StringSet< seqan::CharString >  ids;
                seqan::StringSet< seqan::IupacString > seqs;
                seqan::readRecords( ids, seqs, seqFileIn, config.threads * 5 );
                for ( uint64_t i = 0; i < seqan::length( ids ); ++i )
                {
                    stats.totalSeqsFile += 1;
                    if ( seqan::length( seqs[i] ) < config.kmer_size )
                    { // sequence too small
                        mtx.lock();
                        std::cerr << ids[i] << " has sequence smaller than k-mer size" << std::endl;
                        mtx.unlock();
                        stats.invalidSeqs += 1;
                        continue;
                    }
                    std::string cid   = seqan::toCString( ids[i] );
                    std::string seqid = cid.substr( 0, cid.find( ' ' ) );
                    if ( seq_bin.count( seqid ) == 0 )
                    {
                        mtx.lock();
                        std::cerr << seqid << " not defined on seqid-bin file" << std::endl;
                        mtx.unlock();
                        stats.invalidSeqs += 1;
                        continue;
                    }
                    stats.sumSeqLen += seqan::length( seqs[i] );
                    q_Seqs.push( detail::Seqs{ seqid, seqs[i] } );
                }
            }
            seqan::close( seqFileIn );
        }
        finished = true;
        timeLoadSeq.end();
    } ) );
    for ( auto&& task : tasks )
    {
        task.get();
    }
    timeBuild.end();
    //////////////////////////////


    // Store filter
    timeSaveFilter.start();
    seqan::store( filter, seqan::toCString( config.output_filter_file ) );
    timeSaveFilter.end();
    //////////////////////////////

    timeGanon.end();

    std::cerr << std::endl;
    detail::print_time( config, timeGanon, timeLoadFiles, timeLoadSeq, timeLoadFiles, timeSaveFilter );
    detail::print_stats( stats, config, timeBuild );

    return true;
}
