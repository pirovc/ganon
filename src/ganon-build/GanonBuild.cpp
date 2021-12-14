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
    uint64_t totalBinsBinId;
    uint64_t totalSeqsFile;
    uint64_t totalBinsFile;
    uint64_t invalidSeqs;
    uint64_t newBins;
    double   filterSizeMB;
    double   filterFalsePositive;
};

typedef robin_hood::unordered_map< std::string, std::vector< FragmentBin > > TSeqBin;
typedef std::set< uint64_t >                                                 TBinIds;
typedef std::map< uint64_t, uint64_t >                                       TBinLen;
typedef seqan3::interleaved_bloom_filter<>                                   TFilter;

inline std::string get_seqid( std::string header )
{
    return header.substr( 0, header.find( ' ' ) );
}

void parse_seqid_bin( const std::string& seqid_bin_file, TSeqBin& seq_bin, TBinLen& bin_len )
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
        seq_bin[fields[0]].push_back( FragmentBin{ std::stoul( fields[1] ), std::stoul( fields[2] ), binid } );
        bin_len[binid] += std::stoul( fields[2] ) - std::stoul( fields[1] ) + 1;
    }
}

void parse_sequences( auto& reference_files, TSeqBin& seq_bin, TBinLen& bin_len, uint64_t parse_sequences )
{
    uint64_t binid = parse_sequences;
    for ( auto const& reference_file : reference_files )
    {
        seqan3::sequence_file_input fin{ reference_file };

        for ( auto& [seq, header, qual] : fin )
        {
            // Header id goes up-to first empty space
            const auto seqid = get_seqid( header );
            seq_bin[seqid].push_back( FragmentBin{ 1, seq.size(), binid } );
            bin_len[binid] += seq.size();
        }
        binid++;
    }
}

void parse_refs( SafeQueue< detail::Seqs >& queue_refs,
                 detail::TSeqBin&           seq_bin,
                 std::mutex&                mtx,
                 Stats&                     stats,
                 GanonBuild::Config&        config )
{

    uint16_t wk_size = (config.window_size > 0) ? config.window_size : config.kmer_size;

    for ( auto const& reference_file : config.reference_files )
    {
        // Open file (type define by extension)
        seqan3::sequence_file_input fin{ reference_file };

        // read in chuncks of config.n_refs
        for ( auto&& records : fin | ranges::views::chunk( config.n_refs ) )
        {
            for ( auto& [seq, header, qual] : records )
            {
                const auto seqid = get_seqid( header );
                stats.totalSeqsFile += 1;
                if ( seq.size() < wk_size )
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

                stats.sumSeqLen += seq.size();
                queue_refs.push( detail::Seqs{ std::move( seqid ), std::move( seq ), std::move( fb ) } );
            }
        }
    }
    queue_refs.notify_push_over();
}

template < class Thashes >
uint64_t count_hashes( auto& reference_files, TSeqBin& seq_bin, TBinLen& bin_len, Thashes& hashes_view )
{
    // map instead of vector to easily erase
    robin_hood::unordered_map< uint64_t, robin_hood::unordered_set< uint64_t > > hashes;

    // store parsed lenght. Once whole bin is calculated, remove hashes from memory
    robin_hood::unordered_map< uint64_t, uint64_t > bin_parsed_length;

    uint64_t finished   = 0;
    uint64_t max_hashes = 0;
    uint64_t bin_count  = bin_len.size();
    for ( auto const& reference_file : reference_files )
    {
        seqan3::sequence_file_input fin{ reference_file };

        for ( auto& [seq, header, qual] : fin )
        {
            const auto seqid = get_seqid( header );
            if ( seq_bin.count( seqid ) > 0 )
            {
                for ( const auto& [fragstart, fragend, binid] : seq_bin.at( seqid ) )
                {
                    if ( bin_parsed_length[binid] == bin_len[binid] )
                        continue;

                    const auto mh =
                        seq | seqan3::views::slice( fragstart - 1, fragend ) | hashes_view | std::views::common;
                    hashes[binid].insert( mh.begin(), mh.end() );
                    bin_parsed_length[binid] += ( fragend - fragstart + 1 );

                    if ( bin_parsed_length[binid] == bin_len[binid] )
                    {
                        // if has more hashes than max, keep number
                        if ( hashes[binid].size() > max_hashes )
                            max_hashes = hashes[binid].size();
                        hashes.erase( binid );
                        finished++;
                    }
                }
                if ( finished == bin_count )
                    break;
            }
        }
    }
    // parsed lenght did not achieve max len defined on seqid_bin
    // possible missing sequences or wrong info on seqid_bin
    // count parsed existing hashes
    for ( const auto& [binid, h] : hashes )
    {
        if ( h.size() > max_hashes )
            max_hashes = h.size();
    }

    return max_hashes;
}

inline uint64_t get_optimal_bins( uint64_t nbins )
{
    return std::ceil( nbins / 64.0 ) * 64;
}

inline uint64_t get_bin_size( double false_positive, uint16_t hash_functions, uint64_t max_hashes )
{
    return std::ceil( max_hashes
                      * ( -hash_functions / std::log( 1 - std::exp( std::log( false_positive ) / hash_functions ) ) ) );
}

inline double get_fp( uint64_t bin_size, uint16_t hash_functions, uint64_t max_hashes )
{
    // std::pow( ( 1 - ( std::pow( ( 1 - ( 1 / static_cast< double >( bin_size ) ) ), ( hash_functions * max_hashes ) )
    // ) ), hash_functions );
    return std::pow( 1 - std::exp( -hash_functions / ( bin_size / static_cast< double >( max_hashes ) ) ),
                     hash_functions );
}

void print_ibf_stats(
    Stats& stats, uint64_t max_hashes, uint64_t bin_size, uint64_t bin_count, uint64_t optimal_bin_count )
{
    std::cerr << "IBF stats: " << std::endl;
    std::cerr << "max hashes: " << max_hashes << std::endl;
    std::cerr << "bin count/optimal: " << bin_count << "/" << optimal_bin_count << std::endl;
    std::cerr << "bin size bits: " << bin_size << std::endl;
    std::cerr << "filter size bits/MB: " << ( optimal_bin_count * bin_size ) << "/" << ( stats.filterSizeMB )
              << std::endl;
    std::cerr << "filter max false positive: " << stats.filterFalsePositive << std::endl;
    std::cerr << std::endl;
}

uint64_t get_max_hashes( GanonBuild::Config& config, detail::TSeqBin& seq_bin, TBinLen& bin_len )
{

    uint64_t max_hashes = 0;
    if ( config.count_hashes )
    {
        // count hashes for minimizers or kmers
        if ( config.window_size > 0 )
        {
            auto minimiser_hash = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                                 seqan3::window_size{ config.window_size },
                                                                 seqan3::seed{ 0 } );
            max_hashes          = count_hashes( config.reference_files, seq_bin, bin_len, minimiser_hash );
        }
        else
        {
            auto kmer_hash = seqan3::views::kmer_hash( seqan3::ungapped{ config.kmer_size } );
            max_hashes     = count_hashes( config.reference_files, seq_bin, bin_len, kmer_hash );
        }
    }
    else
    {
        // calculate max hashes based on biggest bin lenght
        auto max_bin_len =
            std::max_element( bin_len.begin(),
                              bin_len.end(),
                              []( const std::pair< uint64_t, uint64_t >& p1,
                                  const std::pair< uint64_t, uint64_t >& p2 ) { return p1.second < p2.second; } );
        max_hashes = max_bin_len->second - config.kmer_size + 1;
    }

    return max_hashes;
}

TFilter create_filter( GanonBuild::Config& config, Stats& stats, uint64_t bin_count, uint64_t max_hashes )
{
    // calculate size/false positive
    double   false_positive    = 0;
    uint64_t bin_size          = 0;
    uint64_t optimal_bin_count = get_optimal_bins( bin_count );
    // Calculate bin size based on filter size (1MB = 8388608bits)
    if ( config.false_positive )
    {
        false_positive = config.false_positive;
        bin_size       = get_bin_size( false_positive, config.hash_functions, max_hashes );
    }
    else
    {
        if ( config.filter_size_mb > 0 )
            bin_size = ( config.filter_size_mb / static_cast< double >( optimal_bin_count ) ) * 8388608u;
        else
            bin_size = config.bin_size_bits;
        false_positive = get_fp( bin_size, config.hash_functions, max_hashes );
    }

    stats.filterSizeMB        = ( optimal_bin_count * bin_size ) / static_cast< double >( 8388608u );
    stats.filterFalsePositive = false_positive;

    if ( config.verbose )
        detail::print_ibf_stats( stats, max_hashes, bin_size, bin_count, optimal_bin_count );

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

void clear_filter( TFilter& filter, const TBinLen& bin_len )
{
    // Reset bins if complete set of sequences is provided (re-create updated bins)
    std::vector< seqan3::bin_index > updated_bins;
    // For all binids in the file provided, only clean bins for the old bins
    // new bins are already cleared by default when created
    for ( auto const& [binid, len] : bin_len )
    {
        if ( binid >= filter.bin_count() - 1 )
        {
            break;
        }
        updated_bins.emplace_back( seqan3::bin_index{ binid } );
    }
    filter.clear( updated_bins );
}

void increase_filter( Config& config, Stats& stats, TFilter& filter, uint64_t bin_count, uint64_t max_hashes )
{
    // Calculate false positive
    uint64_t bin_size          = filter.bin_size();
    uint64_t optimal_bin_count = get_optimal_bins( bin_count );
    double   false_positive    = get_fp( bin_size, config.hash_functions, max_hashes );

    stats.filterSizeMB        = ( optimal_bin_count * bin_size ) / static_cast< double >( 8388608u );
    stats.filterFalsePositive = false_positive;

    if ( config.verbose )
        detail::print_ibf_stats( stats, max_hashes, bin_size, bin_count, optimal_bin_count );

    // If new bins were added
    uint64_t old_total_bins = filter.bin_count();
    stats.newBins           = bin_count - old_total_bins;

    if ( bin_count > filter.bin_count() )
    {
        // just resize if number of bins is bigger than amount on IBF
        // when updating an IBF with empty bins or removing the last bins, this will not be true
        // if new bins are smaller (less bins, sequences removed) IBF still keep all bins but empty
        filter.increase_bin_number_to( seqan3::bin_count{ bin_count } );
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
                    seqan3::debug_stream << hash << ",";
                    filter.emplace( hash, seqan3::bin_index{ binid } );
                }
                seqan3::debug_stream << '\n';
            }
        }
        else
        {
            break;
        }
    }
}

void print_time( const StopClock& timeGanon,
                 const StopClock& timeLoadFiles,
                 const StopClock& timeCountHashes,
                 const StopClock& timeLoadSeq,
                 const StopClock& timeBuild,
                 const StopClock& timeLoadFilter,
                 const StopClock& timeSaveFilter )
{
    using ::operator<<;
    std::cerr << "Loading files     start: " << timeLoadFiles.begin() << std::endl;
    std::cerr << "                    end: " << timeLoadFiles.end() << std::endl;
    std::cerr << "                elapsed: " << timeLoadFiles.elapsed() << std::endl;
    std::cerr << "Counting hashes   start: " << timeCountHashes.begin() << std::endl;
    std::cerr << "                    end: " << timeCountHashes.end() << std::endl;
    std::cerr << "                elapsed: " << timeCountHashes.elapsed() << std::endl;
    std::cerr << "Loading sequences start: " << timeLoadSeq.begin() << std::endl;
    std::cerr << "                    end: " << timeLoadSeq.end() << std::endl;
    std::cerr << "                elapsed: " << timeLoadSeq.elapsed() << std::endl;
    std::cerr << "Loading filter    start: " << timeLoadFilter.begin() << std::endl;
    std::cerr << "                    end: " << timeLoadFilter.end() << std::endl;
    std::cerr << "                elapsed: " << timeLoadFilter.elapsed() << std::endl;
    std::cerr << "Building filter   start: " << timeBuild.begin() << std::endl;
    std::cerr << "                    end: " << timeBuild.end() << std::endl;
    std::cerr << "                elapsed: " << timeBuild.elapsed() << std::endl;
    std::cerr << "Saving filer      start: " << timeSaveFilter.begin() << std::endl;
    std::cerr << "                    end: " << timeSaveFilter.end() << std::endl;
    std::cerr << "                elapsed: " << timeSaveFilter.elapsed() << std::endl;
    std::cerr << "ganon-build       start: " << timeGanon.begin() << std::endl;
    std::cerr << "                    end: " << timeGanon.end() << std::endl;
    std::cerr << "                elapsed: " << timeGanon.elapsed() << std::endl;
    std::cerr << std::endl;
}

void print_stats( Stats& stats, const StopClock& timeGanon )
{
    double   elapsed   = timeGanon.elapsed();
    uint64_t validSeqs = stats.totalSeqsFile - stats.invalidSeqs;
    std::cerr << "ganon-build processed " << validSeqs << " sequences (" << stats.sumSeqLen / 1000000.0 << " Mbp) in "
              << elapsed << " seconds (" << ( validSeqs / 1000.0 ) / ( elapsed / 60.0 ) << " Kseq/m, "
              << ( stats.sumSeqLen / 1000000.0 ) / ( elapsed / 60.0 ) << " Mbp/m)" << std::endl;

    if ( stats.newBins > 0 )
        std::cerr << " - " << stats.newBins << " bins were added to the IBF" << std::endl;
    if ( stats.invalidSeqs > 0 )
        std::cerr << " - " << stats.invalidSeqs << " invalid sequences were skipped" << std::endl;

    std::cerr << " - " << validSeqs << " sequences / " << stats.totalBinsFile << " bins written to the IBF ("
              << std::setprecision( 2 ) << std::fixed << stats.filterSizeMB << "MB)" << std::endl;
    std::cerr << " - " << std::setprecision( 3 ) << std::fixed << stats.filterFalsePositive
              << " max. false positive (only added sequences)" << std::endl;
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
    StopClock       timeCountHashes;
    StopClock       timeLoadFilter;
    StopClock       timeLoadSeq;
    StopClock       timeBuild;
    StopClock       timeSaveFilter;
    std::mutex      mtx;
    detail::Stats   stats;
    detail::TFilter filter;
    detail::TSeqBin seq_bin; // Map with seqid_bin file info
    detail::TBinLen bin_len; // Map with binids used and their total lenght (bp)

    // SafeQueue with reference sequences
    // config.n_batches*config.n_refs = max. amount of references in memory
    SafeQueue< detail::Seqs > queue_refs( config.n_batches );

    // load filter if provided
    timeLoadFilter.start();
    if ( !config.update_filter_file.empty() )
    {
        filter = detail::load_filter( config.update_filter_file );
    }
    timeLoadFilter.stop();

    // load seq_bin and bin_len from --seqid-bin file or create bins automatically
    timeLoadFiles.start();
    if ( !config.seqid_bin_file.empty() )
    {
        detail::parse_seqid_bin( config.seqid_bin_file, seq_bin, bin_len );
    }
    else
    {
        // Iterate over sequences to get their bins and total bin lenghts
        // if updating, create additional bins starting from the last
        uint64_t start_bin = !config.update_filter_file.empty() ? filter.bin_count() : 0;
        detail::parse_sequences( config.reference_files, seq_bin, bin_len, start_bin );
    }
    timeLoadFiles.stop();

    timeCountHashes.start();
    // last key, last bin (std::map)
    uint64_t last_bin = bin_len.rbegin()->first + 1;
    //  based on counted hashes or length of the largest bin
    uint64_t max_hashes = get_max_hashes( config, seq_bin, bin_len );
    timeCountHashes.stop();

    // If updating: extend and clear filter
    timeLoadFilter.start();
    if ( config.update_filter_file.empty() )
    {
        // create new filter:
        filter = detail::create_filter( config, stats, last_bin, max_hashes );
    }
    else
    {
        detail::increase_filter( config, stats, filter, last_bin, max_hashes );
        if ( config.update_complete )
        {
            detail::clear_filter( filter, bin_len );
        }
    }
    stats.totalBinsFile  = filter.bin_count();
    stats.totalBinsBinId = last_bin;
    timeLoadFilter.stop();

    // Start thread for reading the input reference files
    timeLoadSeq.start();
    std::future< void > read_task = std::async( std::launch::async,
                                                detail::parse_refs,
                                                std::ref( queue_refs ),
                                                std::ref( seq_bin ),
                                                std::ref( mtx ),
                                                std::ref( stats ),
                                                std::ref( config ) );

    // Start threads to build filter
    timeBuild.start();
    std::vector< std::future< void > > tasks;
    if ( config.window_size > 0 )
    {
        for ( uint16_t taskNo = 0; taskNo < config.threads_build; ++taskNo )
        {
            auto minimiser_hash = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                                 seqan3::window_size{ config.window_size }
                                                                 ); //, seqan3::seed{ 0 } );
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
        if ( config.verbose )
        {
            detail::print_time(
                timeGanon, timeLoadFiles, timeCountHashes, timeLoadSeq, timeLoadFiles, timeLoadFilter, timeSaveFilter );
        }
        detail::print_stats( stats, timeGanon );
    }
    return true;
}

} // namespace GanonBuild
