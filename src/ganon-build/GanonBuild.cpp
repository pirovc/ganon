#include "GanonBuild.hpp"

#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <cereal/archives/binary.hpp>


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
    std::string                 seqid;
    std::vector< seqan3::dna5 > seq;
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

typedef seqan3::interleaved_bloom_filter<> Tfilter;

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

void store_filter( Tfilter const& filter, std::string const& output_filter_file )
{
    std::ofstream               os( output_filter_file, std::ios::binary ); // Where output should be stored.
    cereal::BinaryOutputArchive archive( os ); // Create an output archive from the output stream.
    archive( filter );                         // Store data.
}


void load_filter( Tfilter& filter, GanonBuild::Config& config, const std::set< uint64_t >& bin_ids, Stats& stats )
{

    // load filter
    std::ifstream              is( config.update_filter_file, std::ios::binary ); // Where input can be found.
    cereal::BinaryInputArchive archive( is ); // Create an input archive from the input stream.
    archive( filter );                        // Load data.

    stats.totalBinsFile = filter.bin_count();

    // last element (set is ordered) plus one
    uint32_t number_new_bins = *bin_ids.rbegin() + 1;
    if ( number_new_bins > stats.totalBinsFile )
    {
        // just resize if number of bins is bigger than amount on IBF
        // when updating an IBF with empty bins or removing the last bins, this will not be true
        filter.increase_bin_number_to( seqan3::bin_count{ number_new_bins } );
        stats.newBins = number_new_bins - stats.totalBinsFile;
    } // if new bins are smaller (less bins, sequences removed) IBF still keep all bins but empty

    // Reset bins if complete set of sequences is provided (re-create updated bins)
    if ( config.update_complete )
    {
        std::vector< seqan3::bin_index > updated_bins;
        // For all binids in the file provided, only clean bins for the old bins
        // new bins are already cleared
        for ( auto const& binid : bin_ids )
        {
            if ( binid >= stats.totalBinsFile - 1 )
            {
                break;
            }
            updated_bins.emplace_back( seqan3::bin_index{ binid } );
        }
        filter.clear( updated_bins );
    }
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
            // Open file (type define by extension)
            seqan3::sequence_file_input fin{ reference_file };

            // read in chuncks of config.n_refs
            for ( auto&& records : fin | ranges::views::chunk( config.n_refs ) )
            {
                for ( auto& [seq, id, qual] : records )
                {
                    stats.totalSeqsFile += 1;
                    if ( seq.size() < config.kmer_size )
                    {
                        if ( config.verbose )
                        {
                            std::scoped_lock lock( mtx );
                            std::cerr << "WARNING: sequence smaller than k-mer size"
                                      << " [" << id << "]" << std::endl;
                        }
                        stats.invalidSeqs += 1;
                        continue;
                    }

                    // Header id goes up-to first empty space
                    std::string seqid = id.substr( 0, id.find( ' ' ) );

                    if ( seq_bin.count( seqid ) == 0 )
                    {
                        if ( config.verbose )
                        {
                            std::scoped_lock lock( mtx );
                            std::cerr << "WARNING: sequence not defined on seqid-bin-file"
                                      << " [" << seqid << "]" << std::endl;
                        }
                        stats.invalidSeqs += 1;
                        continue;
                    }
                    stats.sumSeqLen += seq.size();
                    queue_refs.push( detail::Seqs{ std::move( seqid ), std::move( seq ) } );
                }
            }
        }
        queue_refs.notify_push_over();
    } ) );

    // load new or given filter
    timeLoadFilter.start();
    detail::Tfilter filter;
    if ( !config.update_filter_file.empty() )
    {
        load_filter( filter, config, bin_ids, stats );
    }
    else
    {
        uint64_t bsize;
        // Calculate bin size based on filter size (1MB = 8388608bits)
        if ( config.bin_size_bits == 0 && config.filter_size_mb > 0 )
        {
            uint64_t optimal_bins = ( std::floor( stats.totalBinsBinId / 64 ) + 1 ) * 64;
            bsize                 = ( config.filter_size_mb / static_cast< float >( optimal_bins ) ) * 8388608u;
            std::cerr << "optimal_bins:" << optimal_bins << std::endl;
        }
        else
        {
            bsize = config.bin_size_bits;
        }

        // New filter
        filter = detail::Tfilter{ seqan3::bin_count{ stats.totalBinsBinId },
                                  seqan3::bin_size{ bsize },
                                  seqan3::hash_function_count{ config.hash_functions } };
    }
    stats.totalBinsFile = filter.bin_count();

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
                    for ( uint64_t i = 0; i < seq_bin.at( val.seqid ).size(); i++ )
                    {
                        auto [fragstart, fragend, binid] = seq_bin.at( val.seqid )[i];
                        // Fragment sequence and generate hashes for k-mers
                        auto hashes =
                            val.seq | seqan3::views::slice( fragstart - 1, fragend )
                            | seqan3::views::kmer_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } } );
                        for ( auto const& hash : hashes )
                        {
                            filter.emplace( hash, seqan3::bin_index{ binid } );
                        }
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
    detail::store_filter( filter, config.output_filter_file );
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
