#include "GanonBuild.hpp"

#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <cereal/archives/binary.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

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

struct FragmentBin
{
    uint64_t start;
    uint64_t end;
    uint64_t bin_id;
};

struct Seqs
{
    std::string                 seqid;
    std::vector< seqan3::dna5 > seq;
    std::vector< FragmentBin >  fragbin;
};

struct Stats
{
    Stats()
    : sumSeqLen{ 0 }
    , totalBinsBinId{ 0 }
    , totalSeqsFile{ 0 }
    , totalBinsFile{ 0 }
    , invalidSeqs{ 0 }
    , newBins{ 0 }
    {
    }

    uint64_t sumSeqLen;
    uint32_t totalBinsBinId;
    uint64_t totalSeqsFile;
    uint32_t totalBinsFile;
    uint64_t invalidSeqs;
    uint32_t newBins;
};

typedef std::map< std::string, std::vector< FragmentBin > > TSeqBin;
typedef std::set< uint64_t >                                TBinIds;
typedef seqan3::interleaved_bloom_filter<>                  TFilter;

void parse_seqid_bin( const std::string& seqid_bin_file, TSeqBin& seq_bin, TBinIds& bin_ids )
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
        // save list of unique bins used
        bin_ids.insert( binid );
    }
}

void parse_refs( SafeQueue< detail::Seqs >& queue_refs,
                 detail::TSeqBin&           seq_bin,
                 std::mutex&                mtx,
                 Stats&                     stats,
                 GanonBuild::Config&        config,
                 uint64_t                   nbin_aux )
{

    for ( auto const& reference_file : config.reference_files )
    {
        // Open file (type define by extension)
        seqan3::sequence_file_input fin{ reference_file };

        // read in chuncks of config.n_refs
        for ( auto&& records : fin | ranges::views::chunk( config.n_refs ) )
        {
            for ( auto& [seq, id, qual] : records )
            {
                // Header id goes up-to first empty space
                std::string seqid = id.substr( 0, id.find( ' ' ) );

                stats.totalSeqsFile += 1;
                if ( seq.size() < config.kmer_size )
                {
                    if ( config.verbose )
                    {
                        std::scoped_lock lock( mtx );
                        std::cerr << "WARNING: skipping sequence smaller than k-mer size"
                                  << " [" << seqid << "]" << std::endl;
                    }
                    stats.invalidSeqs += 1;
                    continue;
                }

                std::vector< detail::FragmentBin > fb;
                if ( seq_bin.empty() )
                {
                    // insert whole sequence coordinates into the bin for the reference file
                    // nbin_aux is set with the bin number to start adding sequences
                    // 0 for build and first bin for update
                    fb.push_back( detail::FragmentBin{ 1, seq.size(), nbin_aux } );
                }
                else
                {
                    // get fragments in the seq_bin
                    if ( seq_bin.count( seqid ) > 0 )
                    {
                        fb = seq_bin.at( seqid );
                    }
                    else
                    {
                        // seqid not found in the seq_bin
                        if ( config.verbose )
                        {
                            std::scoped_lock lock( mtx );
                            std::cerr << "WARNING: skipping sequence not defined on seqid-bin-file"
                                      << " [" << seqid << "]" << std::endl;
                        }
                        stats.invalidSeqs += 1;
                        continue;
                    }
                }
                stats.sumSeqLen += seq.size();
                queue_refs.push( detail::Seqs{ std::move( seqid ), std::move( seq ), std::move( fb ) } );
            }
        }
        nbin_aux++;
    }
    queue_refs.notify_push_over();
}

TFilter create_filter( GanonBuild::Config& config, uint32_t bcount )
{
    uint64_t bsize;
    // Calculate bin size based on filter size (1MB = 8388608bits)
    if ( config.bin_size_bits == 0 && config.filter_size_mb > 0 )
    {
        uint64_t optimal_bins = ( std::floor( bcount / 64 ) + 1 ) * 64;
        bsize                 = ( config.filter_size_mb / static_cast< float >( optimal_bins ) ) * 8388608u;
    }
    else
    {
        bsize = config.bin_size_bits;
    }

    return TFilter{ seqan3::bin_count{ bcount },
                    seqan3::bin_size{ bsize },
                    seqan3::hash_function_count{ config.hash_functions } };
}

TFilter load_filter( std::string const& input_filter_file )
{

    TFilter                    filter;
    std::ifstream              is( input_filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( filter );
    return filter;
}

void clear_filter( TFilter& filter, const TBinIds& bin_ids )
{
    // Reset bins if complete set of sequences is provided (re-create updated bins)
    std::vector< seqan3::bin_index > updated_bins;
    // For all binids in the file provided, only clean bins for the old bins
    // new bins are already cleared by default when created
    for ( auto const& binid : bin_ids )
    {
        if ( binid >= filter.bin_count() - 1 )
        {
            break;
        }
        updated_bins.emplace_back( seqan3::bin_index{ binid } );
    }
    filter.clear( updated_bins );
}

void increase_filter( TFilter& filter, uint32_t new_total_bins )
{
    // If new bins were added
    if ( new_total_bins > filter.bin_count() )
    {
        // just resize if number of bins is bigger than amount on IBF
        // when updating an IBF with empty bins or removing the last bins, this will not be true
        // if new bins are smaller (less bins, sequences removed) IBF still keep all bins but empty
        filter.increase_bin_number_to( seqan3::bin_count{ new_total_bins } );
    }
}

void save_filter( TFilter const& filter, std::string const& output_filter_file )
{
    std::ofstream               os( output_filter_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );
    archive( filter );
}

void build( TFilter& filter, SafeQueue< detail::Seqs >& queue_refs, GanonBuild::Config const& config )
{
    auto hash_adaptor = seqan3::views::kmer_hash( seqan3::ungapped{ config.kmer_size } );
    while ( true )
    {
        detail::Seqs val = queue_refs.pop();
        if ( val.seqid != "" )
        {
            for ( uint64_t i = 0; i < val.fragbin.size(); i++ )
            {
                // Fragment sequences
                auto [fragstart, fragend, binid] = val.fragbin[i];
                for ( auto&& hash : val.seq | seqan3::views::slice( fragstart - 1, fragend ) | hash_adaptor )
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
    std::cerr << " - " << validSeqs << " sequences in " << stats.totalBinsFile << " bins written to the IBF"
              << std::endl;
}

} // namespace detail

bool run( Config config )
{

    // Validate configuration input
    if ( !config.validate() )
        return false;

    // Print config
    if ( config.verbose )
        std::cerr << config;

    // Start time count
    StopClock timeGanon;
    timeGanon.start();

    // Initialize variables
    StopClock       timeLoadFiles;
    StopClock       timeLoadFilter;
    StopClock       timeLoadSeq;
    StopClock       timeBuild;
    StopClock       timeSaveFilter;
    std::mutex      mtx;
    detail::Stats   stats;
    detail::TFilter filter;
    detail::TSeqBin seq_bin; // Map with seqid_bin file info
    detail::TBinIds bin_ids; // Set with binids used

    // SafeQueue with reference sequences
    // config.n_batches*config.n_refs = max. amount of references in memory
    SafeQueue< detail::Seqs > queue_refs( config.n_batches );

    // load filter if provided (do it early to get number of bins)
    timeLoadFilter.start();
    if ( !config.update_filter_file.empty() )
    {
        filter = detail::load_filter( config.update_filter_file );
    }
    timeLoadFilter.stop();

    // load seq_bin and bin_ids from --seqid-bin file or create bins automatically
    timeLoadFiles.start();
    if ( !config.seqid_bin_file.empty() )
    {
        detail::parse_seqid_bin( config.seqid_bin_file, seq_bin, bin_ids );
    }
    else
    {
        // seq_bin should remain empty
        int start_bin = 0;
        if ( !config.update_filter_file.empty() )
        {
            start_bin = filter.bin_count(); // if updating, create additional bins starting from the last
        }
        bin_ids = std::views::iota( start_bin, start_bin + int( config.reference_files.size() ) )
                  | seqan3::views::to< detail::TBinIds >; // fill bin_ids with range
    }
    timeLoadFiles.stop();

    // If updating: extend and clear filter
    timeLoadFilter.start();
    if ( !config.update_filter_file.empty() )
    {
        uint32_t old_total_bins = filter.bin_count();
        uint32_t new_total_bins = *bin_ids.rbegin() + 1; // new total is the last bin (set is ordered) + 1
        detail::increase_filter( filter, new_total_bins );
        stats.newBins = new_total_bins - old_total_bins;
        if ( config.update_complete )
        {
            detail::clear_filter( filter, bin_ids );
        }
    }
    else
    {
        // create new filter
        filter = detail::create_filter( config, bin_ids.size() );
    }
    stats.totalBinsFile  = filter.bin_count();
    stats.totalBinsBinId = bin_ids.size();
    timeLoadFilter.stop();

    // Start thread for reading the input reference files
    timeLoadSeq.start();
    // nbin_aux loaded with first bin in case of automatic bin generation
    uint32_t            nbin_aux  = *bin_ids.begin();
    std::future< void > read_task = std::async( std::launch::async,
                                                detail::parse_refs,
                                                std::ref( queue_refs ),
                                                std::ref( seq_bin ),
                                                std::ref( mtx ),
                                                std::ref( stats ),
                                                std::ref( config ),
                                                nbin_aux );

    // Start threads to build filter
    timeBuild.start();
    std::vector< std::future< void > > tasks;
    for ( uint16_t taskNo = 0; taskNo < config.threads_build; ++taskNo )
    {
        tasks.emplace_back( std::async(
            std::launch::async, detail::build, std::ref( filter ), std::ref( queue_refs ), std::ref( config ) ) );
    }
    // Wait until all references sequences are parsed
    read_task.get();
    timeLoadSeq.stop();

    // Wait until all threads building finish
    for ( auto&& task : tasks )
    {
        task.get();
    }
    timeBuild.stop();

    // Store filter
    timeSaveFilter.start();
    detail::save_filter( filter, config.output_filter_file );
    timeSaveFilter.stop();

    // Stop timer for ganon run
    timeGanon.stop();

    // Print time and statistics
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
