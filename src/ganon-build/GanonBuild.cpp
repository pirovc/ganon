#include "GanonBuild.hpp"

#include <robin_hood.h>

#include <defaults/defaults.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>
#include <utils/adjust_seed.hpp>
#include <utils/dna4_traits.hpp>

#include <cereal/archives/binary.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <cinttypes>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

namespace GanonBuild
{

namespace detail
{

typedef robin_hood::unordered_map< std::string, std::string > TTarget;
typedef robin_hood::unordered_map< std::string, uint64_t >    THashesCount;

typedef robin_hood::unordered_map< uint64_t, std::string >                                   TBinMap;
typedef robin_hood::unordered_map< uint64_t, std::tuple< std::string, uint64_t, uint64_t > > TBinMapHash;

typedef seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed > TIBF;

struct IBFConfig
{
    IBFConfig()
    : n_bins{ 0 }
    , max_hashes_bin{ 0 }
    , hash_functions{ 0 }
    , kmer_size{ 0 }
    , window_size{ 0 }
    , bin_size_bits{ 0 }
    , max_fp{ 0 }
    {
    }

    uint64_t n_bins;
    uint64_t max_hashes_bin;
    uint8_t  hash_functions;
    uint8_t  kmer_size;
    uint16_t window_size;
    uint64_t bin_size_bits;
    double   max_fp;
    double   true_max_fp;
    double   true_avg_fp;
};

inline std::ostream& operator<<( std::ostream& stream, const IBFConfig& ibf_config )
{
    constexpr auto newl{ "\n" };
    constexpr auto separator{ "----------------------------------------------------------------------" };

    stream << "ibf_config:" << newl;
    stream << "n_bins         " << ibf_config.n_bins << newl;
    stream << "max_hashes_bin " << ibf_config.max_hashes_bin << newl;
    stream << "hash_functions " << unsigned( ibf_config.hash_functions ) << newl;
    stream << "kmer_size      " << unsigned( ibf_config.kmer_size ) << newl;
    stream << "window_size    " << ibf_config.window_size << newl;
    stream << "bin_size_bits  " << ibf_config.bin_size_bits << newl;
    stream << "max_fp         " << ibf_config.max_fp << newl;
    stream << "true_max_fp    " << ibf_config.true_max_fp << newl;
    stream << "true_avg_fp    " << ibf_config.true_avg_fp << newl;
    stream << separator << newl;
    return stream;
}

struct InputFileMap
{
    std::string file;
    TTarget     targets;
};


struct Total
{
    uint64_t files             = 0;
    uint64_t invalid_files     = 0;
    uint64_t sequences         = 0;
    uint64_t invalid_sequences = 0;
    uint64_t length_bp         = 0;
};

struct Stats
{
    Total total;
    void  add_totals( std::vector< Total > const& totals )
    {
        // add several totals (from threads) to the stats
        for ( auto const& t : totals )
        {
            total.files += t.files;
            total.invalid_files += t.invalid_files;
            total.sequences += t.sequences;
            total.invalid_sequences += t.invalid_sequences;
            total.length_bp += t.length_bp;
        }
    }
};

inline std::string get_seqid( std::string header )
{
    return header.substr( 0, header.find( ' ' ) );
}

robin_hood::unordered_map< std::string, TTarget > parse_input_file( const std::string& input_file,
                                                                    bool               quiet,
                                                                    Stats&             stats )
{
    robin_hood::unordered_map< std::string, TTarget > input_map;
    std::string                                       line;
    std::ifstream                                     infile( input_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field ); // file [<tab> target <tab> sequence]

        const std::string file = fields[0];
        if ( !std::filesystem::exists( file ) || std::filesystem::file_size( file ) == 0 )
        {
            if ( !quiet )
                std::cerr << "input file not found/empty: " << file << std::endl;
            stats.total.invalid_files++;
            continue;
        }

        if ( fields.size() == 1 )
        {
            // skip repeated files
            if ( !input_map.count( file ) )
            {
                // target is the file itself (filename only wihtout path)
                input_map[file][""] = std::filesystem::path( file ).filename();
            }
        }
        else if ( fields.size() == 2 )
        {
            // provided target in the file
            input_map[file][""] = fields[1];
        }
        else
        {
            // sequence as target
            input_map[file][fields[2]] = fields[1];
        }
        stats.total.files++;
    }

    return input_map;
}

void store_hashes( const std::string                            target,
                   const robin_hood::unordered_set< uint64_t >& hashes,
                   const std::string                            tmp_output_folder )
{

    std::filesystem::path outf{ tmp_output_folder };
    outf += target + ".min";
    std::ofstream outfile{ outf, std::ios::binary };
    for ( auto&& h : hashes )
    {
        outfile.write( reinterpret_cast< const char* >( &h ), sizeof( h ) );
    }
}

void delete_hashes( const std::string target, const std::string tmp_output_folder )
{

    std::filesystem::path outf{ tmp_output_folder };
    outf += target + ".min";
    std::filesystem::remove( outf );
}

void count_hashes( SafeQueue< InputFileMap >& ifm_queue,
                   THashesCount&              hashes_count,
                   const GanonBuild::Config&  config,
                   Total&                     total )
{

    auto minimiser_view = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                         seqan3::window_size{ config.window_size },
                                                         seqan3::seed{ raptor::adjust_seed( config.kmer_size ) } );

    while ( true )
    {
        // Wait here until reads are available or push is over and queue is empty
        InputFileMap ifm = ifm_queue.pop();

        // If batch is empty exit thread
        if ( ifm.file == "" )
            break;

        robin_hood::unordered_map< std::string, robin_hood::unordered_set< uint64_t > > hashes;

        seqan3::sequence_file_input< raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq > > fin{
            ifm.file
        };

        for ( auto& [header, seq] : fin )
        {
            std::string target;
            // file as target
            if ( ifm.targets.count( "" ) )
            {
                target = ifm.targets[""];
            }
            else
            {
                // hashes for each sequence
                const auto seqid = get_seqid( header );
                if ( ifm.targets.count( seqid ) )
                {
                    target = ifm.targets[seqid];
                }
                else
                {
                    total.invalid_sequences++;
                    continue;
                }
            }
            total.sequences++;
            total.length_bp += seq.size();
            const auto mh = seq | minimiser_view | std::views::common;
            hashes[target].insert( mh.begin(), mh.end() );
            hashes_count[target] = hashes[target].size();
            detail::store_hashes( target, hashes[target], config.tmp_output_folder );
        }
    }
}


std::vector< uint64_t > load_hashes( std::string file )
{
    uint64_t                hash;
    std::vector< uint64_t > hashes;
    std::ifstream           infile{ file, std::ios::binary };
    while ( infile.read( reinterpret_cast< char* >( &hash ), sizeof( hash ) ) )
        hashes.push_back( hash );
    return hashes;
}

void save_filter( const GanonBuild::Config& config,
                  const TIBF&               ibf,
                  const IBFConfig&          ibf_config,
                  const THashesCount&       hashes_count,
                  const TBinMap&            bin_map )
{
    std::ofstream               os( config.output_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );

    // archive( defaults::version_string );
    // archive( ibf_config );
    // archive( hashes_count );
    // archive( bin_map );
    archive( ibf );
}

int64_t bf_size( double max_fp, uint64_t n_hashes )
{
    return std::ceil( n_hashes * std::log( max_fp ) / std::log( 1.0 / std::pow( 2.0, std::log( 2 ) ) ) );
}

int64_t bf_size( double max_fp, uint64_t n_hashes, uint8_t hash_functions )
{
    return std::ceil( -static_cast< double >( n_hashes * hash_functions )
                      / std::log( 1 - std::exp( std::log( max_fp ) / hash_functions ) ) );
}

uint8_t hash_functions_from_ratio( uint64_t bin_size_bits, uint64_t n_hashes )
{
    return static_cast< uint8_t >( std::log( 2 ) * ( bin_size_bits / static_cast< double >( n_hashes ) ) );
}


uint64_t number_of_bins( THashesCount const& hashes_count, uint64_t n_hashes )
{
    uint64_t n_bins = 0;
    for ( auto const& [target, count] : hashes_count )
    {
        n_bins += std::ceil( count / static_cast< double >( n_hashes ) );
    }
    return n_bins;
}


double correction_rate( uint64_t max_split_bins, double max_fp, uint8_t hash_functions )
{
    return std::log( 1.0 - std::exp( std::log( 1.0 - std::pow( 1.0 - max_fp, max_split_bins ) ) / hash_functions ) )
           / std::log( 1.0 - std::exp( std::log( max_fp ) / hash_functions ) );
}

inline uint64_t optimal_bins( uint64_t n_bins )
{
    return std::ceil( n_bins / 64.0 ) * 64;
}

inline double false_positive( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t n_hashes )
{
    return std::pow( 1 - std::exp( -hash_functions / ( bin_size_bits / static_cast< double >( n_hashes ) ) ),
                     hash_functions );
}

std::tuple< double, double > true_false_positive( THashesCount const& hashes_count,
                                                  uint64_t            max_hashes_bin,
                                                  uint64_t            bin_size_bits,
                                                  uint8_t             hash_functions )
{

    // Calculate true fp for each target, considering split into several bins (multiple testing)
    double highest_fp = 0;
    double average_fp = 0;
    for ( auto const& [target, count] : hashes_count )
    {
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( max_hashes_bin ) );
        uint64_t n_hashes_bin  = std::ceil( count / static_cast< double >( n_bins_target ) );
        double tfp = 1 - std::pow( 1.0 - false_positive( bin_size_bits, hash_functions, n_hashes_bin ), n_bins_target );
        if ( tfp > highest_fp )
            highest_fp = tfp;
        average_fp += tfp;
    }
    average_fp = average_fp / static_cast< double >( hashes_count.size() );

    return std::make_tuple( highest_fp, average_fp );
}

void optimal_hashes_size( double const        filter_size,
                          IBFConfig&          ibf_config,
                          THashesCount const& hashes_count,
                          uint8_t const       hash_functions,
                          uint8_t const       max_hash_functions )
{


    uint64_t max_hashes = 0;
    for ( auto const& [target, cnt] : hashes_count )
    {
        if ( cnt > max_hashes )
            max_hashes = cnt;
    }

    double min_fp = 1;
    // total + 1 zero index
    for ( size_t n = max_hashes + 1; n > 100; n -= 100 )
    {
        uint64_t n_hashes = n - 1;
        uint64_t n_bins   = number_of_bins( hashes_count, n_hashes );

        int64_t bin_size_bits          = ( filter_size / static_cast< double >( optimal_bins( n_bins ) ) ) * 8388608u;
        uint8_t optimal_hash_functions = hash_functions;
        if ( optimal_hash_functions == 0 )
            optimal_hash_functions = hash_functions_from_ratio( bin_size_bits, n_hashes );

        if ( optimal_hash_functions > max_hash_functions || optimal_hash_functions == 0 )
            optimal_hash_functions = max_hash_functions;

        // using avg n_hashes per split bin
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );
        double   fp             = false_positive( bin_size_bits, optimal_hash_functions, n_hashes );
        double   real_fp        = 1 - std::pow( 1.0 - fp, max_split_bins );

        if ( real_fp < min_fp )
        {
            min_fp                    = real_fp;
            ibf_config.max_hashes_bin = n_hashes;
            ibf_config.bin_size_bits  = bin_size_bits;
            ibf_config.hash_functions = optimal_hash_functions;
            ibf_config.max_fp         = real_fp;
            ibf_config.n_bins         = n_bins;
        }
    }
}

void optimal_hashes_fp( double const        max_fp,
                        IBFConfig&          ibf_config,
                        THashesCount const& hashes_count,
                        uint8_t const       hash_functions,
                        uint8_t const       max_hash_functions )
{
    uint64_t min_filter_size = 0;

    uint64_t max_hashes = 0;
    for ( auto const& [target, cnt] : hashes_count )
    {
        if ( cnt > max_hashes )
            max_hashes = cnt;
    }


    ibf_config.max_fp = max_fp;

    // total + 1 zero index
    for ( size_t n = max_hashes + 1; n > 100; n -= 100 )
    {

        uint64_t n_hashes       = n - 1;
        uint64_t n_bins         = number_of_bins( hashes_count, n_hashes );
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );

        int64_t bin_size_bits;
        uint8_t optimal_hash_functions = hash_functions;
        if ( optimal_hash_functions == 0 )
        {
            bin_size_bits          = bf_size( max_fp, n_hashes );
            optimal_hash_functions = hash_functions_from_ratio( bin_size_bits, n_hashes );
        }
        else
        {
            bin_size_bits = bf_size( max_fp, n_hashes, optimal_hash_functions );
        }

        if ( optimal_hash_functions > max_hash_functions || optimal_hash_functions == 0 )
            optimal_hash_functions = max_hash_functions;

        // Approximate real false positive based on average n_hashes per split bin on the max target
        // to not overestimate the correction rate if bins are not completely full
        uint64_t avg_n_hashes = std::ceil( max_hashes / static_cast< double >( max_split_bins ) );
        double   approx_fp    = false_positive( bin_size_bits, optimal_hash_functions, avg_n_hashes );
        if ( approx_fp > max_fp )
            approx_fp = max_fp;

        double crate = correction_rate( max_split_bins, approx_fp, optimal_hash_functions );

        uint64_t filter_size_bits = bin_size_bits * optimal_bins( n_bins ) * crate;

        if ( filter_size_bits < min_filter_size || min_filter_size == 0 )
        {
            min_filter_size           = filter_size_bits;
            ibf_config.max_hashes_bin = n_hashes;
            ibf_config.bin_size_bits  = bin_size_bits * crate;
            ibf_config.n_bins         = n_bins;
            ibf_config.hash_functions = optimal_hash_functions;
        }
    }
}

TBinMapHash create_bin_map_hash( IBFConfig const& ibf_config, THashesCount const& hashes_count )
{

    uint64_t    binno = 0;
    TBinMapHash bin_map_hash;
    for ( auto const& [target, count] : hashes_count )
    {
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( ibf_config.max_hashes_bin ) );
        uint64_t n_hashes_bin  = std::ceil( count / static_cast< double >( n_bins_target ) );

        if ( n_hashes_bin > ibf_config.max_hashes_bin )
            n_hashes_bin = ibf_config.max_hashes_bin;

        for ( uint64_t i = 0; i < n_bins_target; ++i )
        {
            uint64_t hashes_idx_st = i * n_hashes_bin;
            uint64_t hashes_idx_en = hashes_idx_st + n_hashes_bin - 1;
            if ( hashes_idx_st >= count )
                hashes_idx_en = count - 1;
            bin_map_hash[binno] = std::make_tuple( target, hashes_idx_st, hashes_idx_en );
            binno++;
        }
    }
    return bin_map_hash;
}

void build( TIBF&                       ibf,
            std::atomic< std::size_t >& bin_batches,
            const uint64_t              max_batch,
            const uint64_t              batch_size,
            const TBinMapHash&          bin_map_hash,
            std::string                 tmp_output_folder )
{
    while ( true )
    {
        // Add to atomic bin_batches and store
        uint64_t batch = bin_batches++;
        if ( batch >= max_batch )
            break;

        // Set and check boundaries of batches
        uint64_t batch_start = batch * batch_size;
        uint64_t batch_end   = batch_start + batch_size - 1;
        if ( batch_end > bin_map_hash.size() - 1 )
            batch_end = bin_map_hash.size() - 1;

        robin_hood::unordered_map< std::string, std::vector< uint64_t > > target_hashes;
        // Insert hashes by index to the ibf
        for ( uint64_t binno = batch_start; binno <= batch_end; binno++ )
        {
            auto [target, hashes_idx_st, hashes_idx_en] = bin_map_hash.at( binno );
            // read files just once and store
            if ( !target_hashes.count( target ) )
            {
                auto file             = tmp_output_folder + target + ".min";
                target_hashes[target] = load_hashes( file );
            }
            for ( uint64_t pos = hashes_idx_st; pos <= hashes_idx_en; pos++ )
            {
                ibf.emplace( target_hashes[target][pos], seqan3::bin_index{ binno } );
            }
        }
    }
}

void print_stats( Stats& stats, const IBFConfig& ibf_config, const StopClock& timeGanon )
{
    double elapsed = timeGanon.elapsed();
    std::cerr << "ganon-build processed " << stats.total.sequences << " sequences / " << stats.total.files << " files ("
              << stats.total.length_bp / 1000000.0 << " Mbp) in " << elapsed << " seconds ("
              << ( stats.total.length_bp / 1000000.0 ) / ( elapsed / 60.0 ) << " Mbp/m)" << std::endl;

    if ( stats.total.invalid_files > 0 )
        std::cerr << " - " << stats.total.invalid_files << " invalid files skipped" << std::endl;
    if ( stats.total.invalid_sequences > 0 )
        std::cerr << " - " << stats.total.invalid_sequences << " invalid sequences skipped" << std::endl;

    std::cerr << std::fixed << std::setprecision( 4 ) << " - max. false positive: " << ibf_config.true_max_fp;
    std::cerr << std::fixed << std::setprecision( 4 ) << " (avg.: " << ibf_config.true_avg_fp << ")" << std::endl;
    std::cerr << std::fixed << std::setprecision( 2 ) << " - filter size: "
              << ( optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits ) / static_cast< double >( 8388608u )
              << "MB" << std::endl;
}

} // namespace detail

bool run( Config config )
{

    // Validate configuration input
    if ( !config.validate() )
        return false;

    // Print verbose arguments
    if ( config.verbose )
        std::cerr << config;

    // Start timer for total build time
    StopClock timeGanon;
    timeGanon.start();

    // Init. Stats and Totals (one for each thread)
    detail::Stats                stats;
    std::vector< detail::Total > totals( config.threads );

    // create IBF configuration and set-up fixedparameters
    detail::IBFConfig ibf_config;
    ibf_config.kmer_size   = config.kmer_size;
    ibf_config.window_size = config.window_size;

    // Store number of hashes for each target
    detail::THashesCount hashes_count;

    // Parse valid input file into a queue of InputFileMap
    SafeQueue< detail::InputFileMap > ifm_queue;
    for ( auto const& [file, targets] : detail::parse_input_file( config.input_file, config.quiet, stats ) )
    {
        // Initialize hashes_count for each target to be able to update it threads mode
        for ( auto const& target : targets )
            hashes_count[target.second] = 0;

        // Add to File and Targets to the SafeQueue
        ifm_queue.push( detail::InputFileMap{ file, std::move( targets ) } );
    }
    ifm_queue.notify_push_over();

    if ( ifm_queue.empty() )
    {
        std::cerr << "No valid input files" << std::endl;
        return false;
    }
    else
    {
        stats.total.files = ifm_queue.size();
    }

    // Create temporary output folder if not existing
    if ( config.tmp_output_folder != "" && !std::filesystem::exists( config.tmp_output_folder ) )
        std::filesystem::create_directory( config.tmp_output_folder );

    std::vector< std::future< void > > tasks_count;
    for ( uint16_t taskn = 0; taskn < config.threads; ++taskn )
    {
        tasks_count.emplace_back( std::async( std::launch::async,
                                              detail::count_hashes,
                                              std::ref( ifm_queue ),
                                              std::ref( hashes_count ),
                                              std::ref( config ),
                                              std::ref( totals[taskn] ) ) );
    }
    for ( auto&& task : tasks_count )
    {
        task.get();
    }

    // Define optimal parameters and fill ibf_config
    if ( config.filter_size > 0 )
    {
        // Optimal max hashes per bin based on filter size (smallest fp)
        detail::optimal_hashes_size(
            config.filter_size, ibf_config, hashes_count, config.hash_functions, config.max_hash_functions );
    }
    else
    {
        // Optimal max hashes per bin based on max_fp (smallest filter size)
        detail::optimal_hashes_fp(
            config.max_fp, ibf_config, hashes_count, config.hash_functions, config.max_hash_functions );
    }

    // Split hashes into optimal size creating technical bins
    // {binno: (target, idx_hashes_start, idx_hashes_end)}
    const detail::TBinMapHash bin_map_hash = create_bin_map_hash( ibf_config, hashes_count );

    // Calculate true fp for the choosen parameters
    // when calculating optimal parameters, assumes full bins for the sake of speed
    // in reality bins are rarely completely fully loaded
    std::tie( ibf_config.true_max_fp, ibf_config.true_avg_fp ) = detail::true_false_positive(
        hashes_count, ibf_config.max_hashes_bin, ibf_config.bin_size_bits, ibf_config.hash_functions );

    // Print verbose arguments for ibf
    if ( config.verbose )
        std::cerr << ibf_config;

    // TODO assert bin_map_hash.size() == ibf_config.n_bins

    std::cerr << bin_map_hash.size() << std::endl;
    std::cerr << ibf_config.n_bins << std::endl;

    // Create filter
    auto ibf = detail::TIBF{ seqan3::bin_count{ ibf_config.n_bins },
                             seqan3::bin_size{ ibf_config.bin_size_bits },
                             seqan3::hash_function_count{ ibf_config.hash_functions } };

    // build ibf reading .min files
    // parallelize for every 64 entries in bin_map_hash, insert to ibf
    uint64_t batch_size = 64;
    uint64_t max_batch  = std::ceil( bin_map_hash.size() / static_cast< double >( batch_size ) );

    // Keep track of batches processed in the threads
    std::atomic< std::size_t >         bin_batches = 0;
    std::vector< std::future< void > > tasks_build;
    for ( uint16_t taskn = 0; taskn < config.threads; ++taskn )
    {
        tasks_build.emplace_back( std::async( std::launch::async, [=, &ibf, &bin_batches, &bin_map_hash]() {
            detail::build( ibf, bin_batches, max_batch, batch_size, bin_map_hash, config.tmp_output_folder );
        } ) );
    }
    for ( auto&& task : tasks_build )
    {
        task.get();
    }

    // delete temporary hash files
    for ( auto const& [target, cnt] : hashes_count )
    {
        detail::delete_hashes( target, config.tmp_output_folder );
    }

    detail::TBinMap bin_map;
    std::ofstream   outfile{ config.output_file + ".map" }; // TEMP test - TODO remove
    // Create a map for each bin {binno: target}
    for ( auto& [binno, vals] : bin_map_hash )
    {
        bin_map[binno] = std::get< 0 >( vals );
        outfile << bin_map[binno] << '\t' << binno << std::endl; // TEMP test - TODO remove
    }
    outfile.close(); // TEMP test - TODO remove

    // write ibf and other infos
    detail::save_filter( config, ibf, ibf_config, hashes_count, bin_map );

    // sum totals from threads
    stats.add_totals( totals );

    // Stop timer for total build time
    timeGanon.stop();

    if ( !config.quiet )
        detail::print_stats( stats, ibf_config, timeGanon );

    return true;
}

} // namespace GanonBuild
