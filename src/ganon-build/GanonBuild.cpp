#include "GanonBuild.hpp"

#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <seqan/binning_directory.h>

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

namespace GanonBuild
{

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

typedef std::map< std::string, std::vector< FragmentBin > > TSeqBin;

typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                 seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
    Tfilter;

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


Tfilter load_filter( GanonBuild::Config& config, const std::set< uint64_t >& bin_ids, Stats& stats )
{

    Tfilter filter;

    // load from disk in case of update
    if ( !config.update_filter_file.empty() )
    {
        // load filter
        seqan::retrieve( filter, seqan::toCString( config.update_filter_file ) );
        // totalBinsFile account for all bins, even empty
        stats.totalBinsFile = seqan::getNumberOfBins( filter );
        config.kmer_size    = seqan::getKmerSize( filter );
        // config.hash_functions = seqan::get...( filter ); // not avail.
        // config.filter_size_bits = seqan::get...( filter ); // not avail.

        // last element (set is ordered) plus one

        uint32_t number_new_bins = *bin_ids.rbegin() + 1;
        if ( number_new_bins > stats.totalBinsFile )
        {
            // just resize if number of bins is bigger than amount on IBF
            // when updating an IBF with empty bins or removing the last bins, this will not be true
            filter.resizeBins( number_new_bins );
            stats.newBins = number_new_bins - stats.totalBinsFile;
        } // if new bins are smaller (less bins, sequences removed) IBF still keep all bins but empty

        // Reset bins if complete set of sequences is provided (re-create updated bins)
        if ( config.update_complete )
        {
            std::vector< uint32_t > updated_bins;
            // For all binids in the file provided, only clean bins for the old bins
            // new bins are already cleared
            for ( auto const& binid : bin_ids )
            {
                if ( binid >= stats.totalBinsFile - 1 )
                {
                    break;
                }
                updated_bins.emplace_back( binid );
                // std::cerr << "Cleared: " << binid << std::endl;
            }
            seqan::clear( filter, updated_bins, config.threads ); // clear modified bins
        }
    }
    else
    {
        filter = Tfilter( stats.totalBinsBinId, config.hash_functions, config.kmer_size, config.filter_size_bits );
        stats.totalBinsFile = seqan::getNumberOfBins( filter );
    }


    return filter;
}

void print_time( const GanonBuild::Config& config,
                 const StopClock&          timeGanon,
                 const StopClock&          timeLoadFiles,
                 const StopClock&          timeLoadSeq,
                 const StopClock&          timeBuild,
                 const StopClock&          timeLoadFilter,
                 const StopClock&          timeSaveFilter )
{
    using ::operator<<;

    std::cerr << "ganon-build       start time: " << timeGanon.begin() << std::endl;
    std::cerr << "Loading files     start time: " << timeLoadFiles.begin() << std::endl;
    std::cerr << "Loading files       end time: " << timeLoadFiles.end() << std::endl;
    std::cerr << "Loading sequences start time: " << timeLoadSeq.begin() << std::endl;
    std::cerr << "Loading filter    start time: " << timeLoadFilter.begin() << std::endl;
    std::cerr << "Loading filter      end time: " << timeLoadFilter.end() << std::endl;
    std::cerr << "Building filter   start time: " << timeBuild.begin() << std::endl;
    std::cerr << "Loading sequences   end time: " << timeLoadSeq.end() << std::endl;
    std::cerr << "Building filter     end time: " << timeBuild.end() << std::endl;
    std::cerr << "Saving filter     start time: " << timeSaveFilter.begin() << std::endl;
    std::cerr << "Saving filter       end time: " << timeSaveFilter.end() << std::endl;
    std::cerr << "ganon-build         end time: " << timeGanon.end() << std::endl;
    std::cerr << std::endl;
    std::cerr << " - loading files: " << timeLoadFiles.elapsed() << std::endl;
    std::cerr << " - loading filter: " << timeLoadFilter.elapsed() << std::endl;
    std::cerr << " - loading sequences (1t): " << timeLoadSeq.elapsed() << std::endl;
    std::cerr << " - building filter (" << config.threads_build << "t): " << timeBuild.elapsed() << std::endl;
    std::cerr << " - saving filter: " << timeSaveFilter.elapsed() << std::endl;
    std::cerr << " - total: " << timeGanon.elapsed() << std::endl;
    std::cerr << std::endl;
}

void print_stats( Stats& stats, const StopClock& timeBuild )
{
    double   elapsed_build = timeBuild.elapsed();
    uint64_t validSeqs     = stats.totalSeqsFile - stats.invalidSeqs;
    std::cerr << "ganon-build processed " << validSeqs << " sequences (" << stats.sumSeqLen / 1000000.0 << " Mbp) in "
              << elapsed_build << " seconds (" << ( validSeqs / 1000.0 ) / ( elapsed_build / 60.0 ) << " Kseq/m, "
              << ( stats.sumSeqLen / 1000000.0 ) / ( elapsed_build / 60.0 ) << " Mbp/m)" << std::endl;
    if ( stats.invalidSeqs > 0 )
        std::cerr << " - " << stats.invalidSeqs << " invalid sequences were skipped" << std::endl;
    if ( stats.newBins > 0 )
        std::cerr << " - " << stats.newBins << " bins were added to the IBF" << std::endl;
    std::cerr << " - " << validSeqs << " sequences in " << stats.totalBinsFile + stats.newBins
              << " bins were written to the IBF" << std::endl;
}

} // namespace detail

bool run( Config config )
{

    if ( !config.validate() )
        return false;

    StopClock timeGanon;
    timeGanon.start();
    StopClock timeLoadFiles;
    StopClock timeLoadFilter;
    StopClock timeLoadSeq;
    StopClock timeBuild;
    StopClock timeSaveFilter;

    if ( config.verbose )
        std::cerr << config;

    //////////////////////////////

    detail::Stats stats;

    timeLoadFiles.start();
    // parse seqid bin
    detail::TSeqBin      seq_bin;
    std::set< uint64_t > bin_ids;
    parse_seqid_bin( config.seqid_bin_file, seq_bin, bin_ids );
    stats.totalSeqsBinId = seq_bin.size();
    stats.totalBinsBinId = bin_ids.size();
    timeLoadFiles.stop();

    //////////////////////////////

    std::mutex                mtx;
    SafeQueue< detail::Seqs > queue_refs(
        config.n_batches ); // config.n_batches*config.n_refs = max. amount of references in memory

    // Start extra thread for reading the input
    timeLoadSeq.start();
    std::future< void > read_task( std::async( std::launch::async, [=, &seq_bin, &queue_refs, &mtx, &stats] {
        for ( auto const& reference_file : config.reference_files )
        {
            seqan::SeqFileIn seqFileIn;
            if ( !seqan::open( seqFileIn, seqan::toCString( reference_file ) ) )
            {
                std::cerr << "ERROR: Unable to open the file: " << reference_file << std::endl;
                continue;
            }
            while ( !seqan::atEnd( seqFileIn ) )
            {

                seqan::StringSet< seqan::CharString > ids;
                seqan::StringSet< seqan::CharString > seqs;

                try
                {
                    seqan::readRecords( ids, seqs, seqFileIn, config.n_refs );
                }
                catch ( seqan::Exception const& e )
                {
                    mtx.lock();
                    std::cerr << "ERROR: Problems parsing the file: " << reference_file << "[" << e.what() << "]"
                              << std::endl;
                    mtx.unlock();
                }

                for ( uint64_t i = 0; i < seqan::length( ids ); ++i )
                {
                    stats.totalSeqsFile += 1;
                    if ( config.verbose && seqan::length( seqs[i] ) < config.kmer_size )
                    { // sequence too small
                        mtx.lock();
                        std::cerr << "WARNING: sequence smaller than k-mer size"
                                  << " [" << ids[i] << "]" << std::endl;
                        mtx.unlock();
                        stats.invalidSeqs += 1;
                        continue;
                    }
                    std::string cid   = seqan::toCString( ids[i] );
                    std::string seqid = cid.substr( 0, cid.find( ' ' ) );
                    if ( config.verbose && seq_bin.count( seqid ) == 0 )
                    {
                        mtx.lock();
                        std::cerr << "WARNING: sequence not defined on seqid-bin-file"
                                  << " [" << seqid << "]" << std::endl;
                        mtx.unlock();
                        stats.invalidSeqs += 1;
                        continue;
                    }
                    stats.sumSeqLen += seqan::length( seqs[i] );
                    queue_refs.push( detail::Seqs{ seqid, seqs[i] } );
                }
            }
            seqan::close( seqFileIn );
        }
        queue_refs.notify_push_over();
    } ) );

    // load new or given filter
    timeLoadFilter.start();
    detail::Tfilter filter = load_filter( config, bin_ids, stats );
    timeLoadFilter.stop();

    // Start execution threads to add kmers
    timeBuild.start();
    std::vector< std::future< void > > tasks;
    for ( uint16_t taskNo = 0; taskNo < config.threads_build; ++taskNo )
    {
        tasks.emplace_back( std::async( std::launch::async, [=, &seq_bin, &filter, &queue_refs] {
            while ( true )
            {
                detail::Seqs val = queue_refs.pop();
                if ( val.seqid != "" )
                {
                    for ( uint64_t i = 0; i < seq_bin[val.seqid].size(); i++ )
                    {
                        auto [fragstart, fragend, binid] = seq_bin[val.seqid][i];
                        // For infixes, we have to provide both the including start and the excluding end position.
                        // fragstart -1 to fix offset
                        // fragend -1+1 to fix offset and not exclude last position
                        seqan::Infix< seqan::Dna5String >::Type fragment = infix( val.seq, fragstart - 1, fragend );
                        seqan::insertKmer( filter, fragment, binid );
                    }
                }
                else
                {
                    break;
                }
            }
        } ) );
    }
    //////////////////////////////

    read_task.get();
    timeLoadSeq.stop();

    for ( auto&& task : tasks )
    {
        task.get();
    }
    timeBuild.stop();
    //////////////////////////////

    // Store filter
    timeSaveFilter.start();
    seqan::store( filter, seqan::toCString( config.output_filter_file ) );
    timeSaveFilter.stop();
    //////////////////////////////

    timeGanon.stop();

    if ( !config.quiet )
    {
        std::cerr << std::endl;
        if ( config.verbose )
        {
            detail::print_time(
                config, timeGanon, timeLoadFiles, timeLoadSeq, timeLoadFiles, timeLoadFilter, timeSaveFilter );
        }
        detail::print_stats( stats, timeBuild );
    }
    return true;
}

} // namespace GanonBuild
