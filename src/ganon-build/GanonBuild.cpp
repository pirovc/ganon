#include "GanonBuild.hpp"

#include <robin_hood.h>

#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <cereal/archives/binary.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

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

void parse_seqid_bin( const std::string&       seqid_bin_file,
                      TSeqBin&                 seq_bin,
                      TBinIds&                 bin_ids,
                      std::vector< uint64_t >& bin_total_length )
{
    std::string   line;
    std::ifstream infile( seqid_bin_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field ); // seqid <tab> seqstart <tab> seqend <tab> binid

        const uint64_t binid = std::stoul( fields[3] );
        if ( binid >= bin_total_length.size() )
        {
            bin_total_length.resize( binid + 1 );
        }

        const uint64_t frag_len = std::stoul( fields[2] ) - std::stoul( fields[1] ) + 1;
        seq_bin[fields[0]].push_back( FragmentBin{ std::stoul( fields[1] ), std::stoul( fields[2] ), binid } );
        bin_total_length[binid] += frag_len;

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
    // min. size necessary to parse reference (either window with minim. or kmer)
    uint8_t base_size = config.window_size > 0 ? config.window_size : config.kmer_size;
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
                if ( seq.size() < base_size )
                {
                    if ( config.verbose )
                    {
                        std::scoped_lock lock( mtx );
                        std::cerr << "WARNING: skipping sequence smaller than k-mer/window size"
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

uint64_t count_hashes( GanonBuild::Config& config, TSeqBin& seq_bin, std::vector< uint64_t >& bin_total_length )
{

    uint64_t bin_count = bin_total_length.size();

    // map instead of vector to easily erase
    robin_hood::unordered_map< uint64_t, robin_hood::unordered_set< uint64_t > > hashes;

    // store parsed lenght. Once whole bin is calculated, remove hashes from memory
    std::vector< uint64_t > bin_parsed_length( bin_count, 0 );

    auto minimiser_hash = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                         seqan3::window_size{ config.window_size },
                                                         seqan3::seed{ 0 } );

    uint64_t finished   = 0;
    uint64_t max_hashes = 0;

    for ( auto const& reference_file : config.reference_files )
    {
        seqan3::sequence_file_input fin{ reference_file };

        for ( auto& [seq, id, qual] : fin )
        {
            // Header id goes up-to first empty space
            const std::string seqid = id.substr( 0, id.find( ' ' ) );

            if ( seq_bin.count( seqid ) > 0 )
            {
                for ( const auto& [fragstart, fragend, binid] : seq_bin.at( seqid ) )
                {
                    if ( bin_parsed_length[binid] == bin_total_length[binid] )
                        continue;

                    const auto mh =
                        seq | seqan3::views::slice( fragstart - 1, fragend ) | minimiser_hash | std::views::common;
                    hashes[binid].insert( mh.begin(), mh.end() );
                    bin_parsed_length[binid] += ( fragend - fragstart + 1 );

                    if ( bin_parsed_length[binid] == bin_total_length[binid] )
                    {
                        if ( hashes[binid].size() > max_hashes )
                            max_hashes = hashes[binid].size();
                        finished++;
                    }
                }
                if ( finished == bin_count )
                    break;
            }
        }
    }
    // parsed lenght did not achieve max len defined on seqid_bin (possible missing sequences)
    // count parsed existing hashes
    for ( const auto& [binid, h] : hashes )
    {
        if ( h.size() > max_hashes )
            max_hashes = h.size();
    }

    return max_hashes;
}

TFilter create_filter( GanonBuild::Config& config, uint32_t bin_count, uint64_t max_hashes )
{

    const auto mb_bits           = 8388608u;
    float      false_positive    = 0;
    uint64_t   bin_size          = 0;
    uint64_t   optimal_bin_count = ( std::floor( bin_count / 64 ) + 1 ) * 64;
    // Calculate bin size based on filter size (1MB = 8388608bits)
    if ( config.false_positive )
    {
        bin_size       = std::ceil( -( 1
                                 / ( pow( ( 1 - pow( config.false_positive, ( 1 / float( config.hash_functions ) ) ) ),
                                          ( 1 / float( config.hash_functions * max_hashes ) ) )
                                     - 1 ) ) );
        false_positive = config.false_positive;
    }
    else
    {
        if ( config.filter_size_mb > 0 )
        {
            bin_size = ( config.filter_size_mb / static_cast< float >( optimal_bin_count ) ) * mb_bits;
        }
        else
        {
            bin_size = config.bin_size_bits;
        }
        false_positive =
            pow( ( 1 - ( pow( ( 1 - ( 1 / float( bin_size ) ) ), ( config.hash_functions * max_hashes ) ) ) ),
                 config.hash_functions );
    }

    if ( config.verbose )
    {
        std::cerr << "max hashes: " << max_hashes << std::endl;
        std::cerr << "bin count/optimal: " << bin_count << "/" << optimal_bin_count << std::endl;
        std::cerr << "bin size bits: " << bin_size << std::endl;
        std::cerr << "filter size bits/MB: " << ( optimal_bin_count * bin_size ) << "/"
                  << ( ( optimal_bin_count * bin_size ) / static_cast< float >( mb_bits ) ) << std::endl;
        std::cerr << "filter max false positive: " << false_positive << std::endl;
    }
    return TFilter{ seqan3::bin_count{ bin_count },
                    seqan3::bin_size{ bin_size },
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

template < class Thashes >
void build( TFilter& filter, SafeQueue< detail::Seqs >& queue_refs, Thashes& hashes_view )
{

    while ( true )
    {
        detail::Seqs val = queue_refs.pop();
        if ( val.seqid != "" )
        {
            for ( uint64_t i = 0; i < val.fragbin.size(); i++ )
            {
                // Fragment sequences
                auto [fragstart, fragend, binid] = val.fragbin[i];
                for ( auto&& hash : val.seq | seqan3::views::slice( fragstart - 1, fragend ) | hashes_view )
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
    StopClock               timeLoadFiles;
    StopClock               timeLoadFilter;
    StopClock               timeLoadSeq;
    StopClock               timeBuild;
    StopClock               timeSaveFilter;
    std::mutex              mtx;
    detail::Stats           stats;
    detail::TFilter         filter;
    detail::TSeqBin         seq_bin;          // Map with seqid_bin file info
    detail::TBinIds         bin_ids;          // Set with binids used
    std::vector< uint64_t > bin_total_length; // Store total lenght (bp) for each bin

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
        detail::parse_seqid_bin( config.seqid_bin_file, seq_bin, bin_ids, bin_total_length );
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
        // Calculae bin_size based on hashes (minimisers) or length of the largest bin (kmers)
        uint64_t max_hashes = 0;
        if ( config.window_size > 0 )
        {
            // count hashes
            max_hashes = count_hashes( config, seq_bin, bin_total_length );
        }
        else
        {
            const auto max_bin_len = *std::max_element( bin_total_length.begin(), bin_total_length.end() );
            max_hashes             = max_bin_len - config.kmer_size + 1;
        }

        // create new filter
        filter = detail::create_filter( config, bin_total_length.size(), max_hashes );
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
    if ( config.window_size > 0 )
    {
        for ( uint16_t taskNo = 0; taskNo < config.threads_build; ++taskNo )
        {
            auto minimiser_hash = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                                 seqan3::window_size{ config.window_size },
                                                                 seqan3::seed{ 0 } );
            tasks.emplace_back( std::async( std::launch::async, [&filter, &queue_refs, &minimiser_hash]() {
                detail::build( filter, queue_refs, minimiser_hash );
            } ) );
        }
    }
    else
    {
        for ( uint16_t taskNo = 0; taskNo < config.threads_build; ++taskNo )
        {
            auto kmer_hash = seqan3::views::kmer_hash( seqan3::ungapped{ config.kmer_size } );
            tasks.emplace_back( std::async( std::launch::async, [&filter, &queue_refs, &kmer_hash]() {
                detail::build( filter, queue_refs, kmer_hash );
            } ) );
        }
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
