#include "GanonBuild.hpp"

#include <robin_hood.h>

#include <defaults/defaults.hpp>
#include <utils/IBFConfig.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>
#include <utils/adjust_seed.hpp>
#include <utils/dna4_traits.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/exception.hpp>
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


typedef seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed > TIBF;

typedef robin_hood::unordered_map< std::string, uint64_t > THashesCount;

typedef robin_hood::unordered_map< uint64_t, std::tuple< std::string, uint64_t, uint64_t > > TBinMapHash;

struct InputFileMap
{
    std::string                target;
    std::vector< std::string > files;
};


struct Total
{
    uint64_t files             = 0;
    uint64_t invalid_files     = 0;
    uint64_t sequences         = 0;
    uint64_t skipped_sequences = 0;
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
            total.skipped_sequences += t.skipped_sequences;
            total.length_bp += t.length_bp;
        }
    }
};

inline std::string get_seqid( std::string header )
{
    /*
     * returns accession version from fasta headers (anything up to first space)
     */
    return header.substr( 0, header.find( ' ' ) );
}

robin_hood::unordered_map< std::string, std::vector< std::string > > parse_input_file( const std::string& input_file,
                                                                                       THashesCount&      hashes_count,
                                                                                       bool               quiet,
                                                                                       Stats&             stats )
{
    /*
     * Funtion to parse input file -> tabular file with the fields: file [<tab> target <tab> seqid]
     * Returns an map with {target: {file: [seqids]}}
     * In case of sequence parsing (provided seqids on third row), seqids is a set of sequence ids
     * In case of file parsing, each file has only one target and seqids is empty
     */
    robin_hood::unordered_map< std::string, std::vector< std::string > > input_map;
    std::string                                                          line;
    robin_hood::unordered_set< std::string >                             files;
    std::ifstream                                                        infile( input_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field ); // file [<tab> target <tab> seqid]

        const std::string file = fields[0];
        files.insert( file );
        if ( !std::filesystem::exists( file ) || std::filesystem::file_size( file ) == 0 )
        {
            if ( !quiet )
                std::cerr << "WARNING: input file not found/empty: " << file << std::endl;
            stats.total.invalid_files++;
            continue;
        }

        if ( fields.size() == 1 )
        {
            // target is the file itself (filename only without path)
            const auto target = std::filesystem::path( file ).filename();
            input_map[target].push_back( file );
            hashes_count[target] = 0;
        }
        else if ( fields.size() == 2 )
        {
            // provided target in the file
            input_map[fields[1]].push_back( file );
            hashes_count[fields[1]] = 0;
        }
    }
    stats.total.files = files.size();

    return input_map;
}


void store_hashes( const std::string                            target,
                   const robin_hood::unordered_set< uint64_t >& hashes,
                   const std::string                            tmp_output_folder )
{
    /*
     * store hashes from set to disk in the specified folder (or current folder ".")
     */
    std::filesystem::path outf{ tmp_output_folder };
    outf += target + ".min";
    std::ofstream outfile{ outf, std::ios::binary | std::ios::app };
    for ( auto&& h : hashes )
    {
        outfile.write( reinterpret_cast< const char* >( &h ), sizeof( h ) );
    }
    outfile.close();
}


std::vector< uint64_t > load_hashes( std::string file )
{
    /*
     * load hashes file from disk and returns them in a vector
     */
    uint64_t                hash;
    std::vector< uint64_t > hashes;
    std::ifstream           infile{ file, std::ios::binary };
    while ( infile.read( reinterpret_cast< char* >( &hash ), sizeof( hash ) ) )
        hashes.push_back( hash );
    return hashes;
}

void delete_hashes( const THashesCount& hashes_count, const std::string tmp_output_folder )
{
    /*
     * delete hashes from disk
     */
    for ( auto const& [target, cnt] : hashes_count )
    {
        std::filesystem::path outf{ tmp_output_folder };
        outf += target + ".min";
        if ( std::filesystem::exists( outf ) )
            std::filesystem::remove( outf );
    }
}

void count_hashes( SafeQueue< InputFileMap >& ifm_queue,
                   THashesCount&              hashes_count,
                   const GanonBuild::Config&  config,
                   Total&                     total )
{

    /*
     * function running in parallel reading from a queue of InputFileMap
     * it counts minimizer hashes for each target and store it to disk
     * it also keep the counts to be further use to define best filter size
     */


    // one view per thread
    auto minimiser_view = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                         seqan3::window_size{ config.window_size },
                                                         seqan3::seed{ raptor::adjust_seed( config.kmer_size ) } );

    while ( true )
    {
        // Wait here until reads are available or push is over and queue is empty
        InputFileMap ifm = ifm_queue.pop();

        // If empty after pop, exit thread
        if ( ifm.target == "" )
            break;

        // For all files of a target
        for ( auto& file : ifm.files )
        {
            try
            {
                // open file
                seqan3::sequence_file_input< raptor::dna4_traits,
                                             seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                    fin{ file };


                // File as target - generate all hashes from file with possible multiple sequences
                // before counting and storing

                robin_hood::unordered_set< uint64_t > hashes;
                for ( auto const& [header, seq] : fin )
                {
                    if ( seq.size() < config.min_length )
                    {
                        total.skipped_sequences++;
                        continue;
                    }

                    total.sequences++;
                    total.length_bp += seq.size();
                    const auto mh = seq | minimiser_view | std::views::common;
                    hashes.insert( mh.begin(), mh.end() );
                }
                hashes_count[ifm.target] += hashes.size();
                detail::store_hashes( ifm.target, hashes, config.tmp_output_folder );
            }
            catch ( seqan3::parse_error const& e )
            {
                std::cerr << "Error parsing file [" << file << "]. " << e.what() << std::endl;
                continue;
            }
        }
    }
}

void save_filter( const GanonBuild::Config& config,
                  const TIBF&               ibf,
                  const IBFConfig&          ibf_config,
                  const THashesCount&       hashes_count,
                  const TBinMapHash&        bin_map_hash )
{
    /*
     * saves structures and IBF to file. structures are coverted to std. lib
     * contents:
     * version:      tuple(major, minor, patch)
     * ibf_config:   IBFConfig struct with values used to build filter
     * hashes_count: vector<tuple(string, int)> [(target, count)]
     * bin_map:      vector<tuple(int, string)> [(binno, target)]
     * ibf:          Interleaved Bloom Filter from seqan3
     */
    std::ofstream               os( config.output_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );

    // Create a map for each bin and store on stdlib structures {binno: target}
    std::vector< std::tuple< uint64_t, std::string > > bin_map;
    for ( auto& [binno, vals] : bin_map_hash )
    {
        bin_map.push_back( { binno, std::get< 0 >( vals ) } );
    }

    // Store hashes_count on stdlib structures {target: hash_count}
    std::vector< std::tuple< std::string, uint64_t > > hashes_count_std;
    for ( auto& [target, count] : hashes_count )
    {
        hashes_count_std.push_back( { target, count } );
    }

    archive( std::make_tuple( defaults::version_major, defaults::version_minor, defaults::version_patch ) );
    archive( ibf_config );
    archive( hashes_count_std );
    archive( bin_map );
    archive( ibf );
}

uint64_t bin_size( double max_fp, uint64_t n_hashes )
{
    /*
     * calculates bin size (bloom filter) based on max. false positive and elements to be inserted
     */
    return std::ceil( ( n_hashes * std::log( max_fp ) ) / std::log( 1.0 / std::pow( 2, std::log( 2 ) ) ) );
}

uint64_t bin_size( double max_fp, uint64_t n_hashes, uint8_t hash_functions )
{
    /*
     * calculates bin size (bloom filter) based on max. false positive, elements to be inserted
     * and number of hash functions
     */
    return std::ceil( n_hashes
                      * ( -hash_functions / std::log( 1 - std::exp( std::log( max_fp ) / hash_functions ) ) ) );
}

uint8_t hash_functions_from_ratio( uint64_t bin_size_bits, uint64_t n_hashes )
{
    /*
     * optimal number of hash functions based on bin size and number of elements
     */
    return static_cast< uint8_t >( std::log( 2 ) * ( bin_size_bits / static_cast< double >( n_hashes ) ) );
}

uint8_t get_optimal_hash_functions( uint64_t bin_size_bits,
                                    uint64_t n_hashes,
                                    uint8_t  hash_functions,
                                    uint8_t  max_hash_functions )
{
    /*
     * Function check if hash should be calculated based on ratio (hash_functions = 0)
     * and if it's on range permitted 1-5
     */

    uint8_t optimal_hash_functions = hash_functions;
    if ( optimal_hash_functions == 0 )
        optimal_hash_functions = hash_functions_from_ratio( bin_size_bits, n_hashes );
    if ( optimal_hash_functions > max_hash_functions || optimal_hash_functions == 0 )
        optimal_hash_functions = max_hash_functions;

    return optimal_hash_functions;
}


uint64_t number_of_bins( THashesCount const& hashes_count, uint64_t n_hashes )
{
    /*
     * number of bins of the IBF given a number of elements (will count bins and split bins)
     */
    uint64_t n_bins = 0;
    for ( auto const& [target, count] : hashes_count )
    {
        n_bins += std::ceil( count / static_cast< double >( n_hashes ) );
    }
    return n_bins;
}


double correction_rate( uint64_t max_split_bins, double max_fp, uint8_t hash_functions, uint64_t n_hashes )
{
    /*
     * calculates the rate that a bin size should increase to accomodate the multiple error
     * problem created by splitting a target into many bins. Based on target splitted amongs
     * the most bins
     */

    double const target_fpr        = 1.0 - std::exp( std::log( 1.0 - max_fp ) / max_split_bins );
    size_t const new_bin_size      = bin_size( target_fpr, n_hashes, hash_functions );
    size_t const original_bin_size = bin_size( max_fp, n_hashes, hash_functions );
    return static_cast< double >( new_bin_size ) / original_bin_size;
}


inline uint64_t optimal_bins( uint64_t n_bins )
{
    /*
     * returns optimal number of bins (multiple of 64) to create the IBF
     */
    return std::ceil( n_bins / 64.0 ) * 64;
}

inline double false_positive( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t n_hashes )
{
    /*
     * calculates the theoretical false positive of a bin (bf) based on parameters
     */
    return std::pow( 1 - std::exp( -hash_functions / ( bin_size_bits / static_cast< double >( n_hashes ) ) ),
                     hash_functions );
}

std::tuple< double, double > true_false_positive( THashesCount const& hashes_count,
                                                  uint64_t            max_hashes_bin,
                                                  uint64_t            bin_size_bits,
                                                  uint8_t             hash_functions )
{
    /*
     * calculates the theoretical false positive (avg and max) of a the IBF
     * based on targets and split bins
     */
    double highest_fp = 0;
    double average_fp = 0;
    // Calculate true fp for each target group, considering split into several bins (multiple testing)
    for ( auto const& [target, count] : hashes_count )
    {
        // Use average number of hashes for each bin to calculate fp
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( max_hashes_bin ) );
        // this can be off by a very small number (rounding ceil on multiple bins)
        uint64_t n_hashes_bin = std::ceil( count / static_cast< double >( n_bins_target ) );

        // false positive for the current target, considering split bins
        double real_fp =
            1.0 - std::pow( 1.0 - false_positive( bin_size_bits, hash_functions, n_hashes_bin ), n_bins_target );

        if ( real_fp > highest_fp )
            highest_fp = real_fp;
        average_fp += real_fp;
    }
    average_fp = average_fp / static_cast< double >( hashes_count.size() );

    return std::make_tuple( highest_fp, average_fp );
}

inline uint64_t get_max_hashes( THashesCount const& hashes_count )
{
    /*
     * return max number of hashes from hashes_count
     */
    uint64_t max_hashes = 0;
    for ( auto const& [target, cnt] : hashes_count )
    {
        if ( cnt > max_hashes )
            max_hashes = cnt;
    }
    return max_hashes;
}

void optimal_hashes( double const        max_fp,
                     double const        filter_size,
                     IBFConfig&          ibf_config,
                     THashesCount const& hashes_count,
                     uint8_t const       hash_functions,
                     uint8_t const       max_hash_functions,
                     std::string const   mode )
{

    /*
     * given a max. false positive or filter_size, iterate over possible capacities for a bin (single bloom filter)
     * and calculate the respective size, considering split bins
     * selects the parameters generating the "best" between smallest filter/fp and n. of bins and fill the
     * ibf_config struct
     */

    // Target with the highest number of minimizers
    uint64_t const max_hashes = get_max_hashes( hashes_count );

    // save minimal values for average
    uint64_t min_filter_size = 0;
    uint64_t min_bins        = 0;
    double   min_fp          = 1;

    // save simulations for average
    struct SimParam
    {
        uint64_t n_hashes;
        uint64_t n_bins;
        uint64_t filter_size_bits;
        double   fp;
    };
    std::vector< SimParam > simulations;

    // simulation on every 100th n. of elements
    size_t iter = 100;
    // check if max_hashes not smaller than iteration
    if ( max_hashes < iter )
        iter = max_hashes;

    // (total + 1) to deal with zero index
    for ( size_t n = max_hashes + 1; n > iter; n -= iter )
    {
        // number of elements to be inserted in a bin
        uint64_t const n_hashes = n - 1;

        // actual number of bins based on targets and elements (not multiple of 64)
        uint64_t const n_bins = number_of_bins( hashes_count, n_hashes );

        int64_t bin_size_bits          = 0;
        uint8_t optimal_hash_functions = 0;
        if ( filter_size )
        {
            // if filter_size is provided, simple calculation
            bin_size_bits = ( filter_size / static_cast< double >( optimal_bins( n_bins ) ) ) * 8388608u;
            optimal_hash_functions =
                get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
        }
        else
        {
            if ( hash_functions == 0 )
            {
                // First define size and than n. hash functions
                bin_size_bits = bin_size( max_fp, n_hashes );
                optimal_hash_functions =
                    get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
            }
            else
            {
                // n. hash functions provided, define best bin size with it
                optimal_hash_functions =
                    get_optimal_hash_functions( bin_size_bits, n_hashes, hash_functions, max_hash_functions );
                bin_size_bits = bin_size( max_fp, n_hashes, optimal_hash_functions );
            }
        }

        // max. times a target is split into several bins number of splits for one target
        uint64_t const max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );

        // Calculate either fp if filter_size is provided or filter_size if max_fp is provided
        double   fp               = 0;
        uint64_t filter_size_bits = 0;
        if ( filter_size )
        {
            // false positive for the current values, considering split bins
            fp =
                1 - std::pow( 1.0 - false_positive( bin_size_bits, optimal_hash_functions, n_hashes ), max_split_bins );

            // save minimal value
            if ( fp < min_fp )
                min_fp = fp;
        }
        else
        {

            // number of elements actually inserted to each bin
            uint64_t const avg_n_hashes = std::ceil( max_hashes / static_cast< double >( max_split_bins ) );
            // Approximate real false positive based on average n_hashes per split bin
            // if not applied, can overestimate the correction rate if bins are not completely "full"
            double approx_fp = false_positive( bin_size_bits, optimal_hash_functions, avg_n_hashes );
            // if approx is higher (precision calculations) set back
            if ( approx_fp > max_fp )
                approx_fp = max_fp;

            // correction rate based on the max. number of splits of a single target
            double const crate = correction_rate( max_split_bins, approx_fp, optimal_hash_functions, n_hashes );
            // apply to the bin size
            bin_size_bits = bin_size_bits * crate;
            // Calculate final filter size
            filter_size_bits = bin_size_bits * optimal_bins( n_bins );

            // values too small or big due to small n_hashes or too high crate break loop
            if ( filter_size_bits == 0 || std::isinf( crate ) )
                break;

            // save minimal value
            if ( filter_size_bits < min_filter_size || min_filter_size == 0 )
                min_filter_size = filter_size_bits;
        }

        // Save simulation values
        simulations.emplace_back( SimParam{ n_hashes, n_bins, filter_size_bits, fp } );

        /*        std::cout << "n_hashes: " << n_hashes << '\t';
                std::cout << "n_bins: " << n_bins << '\t';
                std::cout << "filter_size_bits: " << filter_size_bits << '\t';
                std::cout << "fp: " << fp << '\t';
                std::cout << "hash_functions: " << unsigned( optimal_hash_functions ) << '\n';*/

        // save minimal value
        if ( n_bins < min_bins || min_bins == 0 )
            min_bins = n_bins;
    }

    // Select "optimal" hashes as a harmonic mean of n_bins and filter_size
    // considering their difference to possible minimal values
    // if special mode is selected, mean is deviated by a factor (smaller is better for ratios, so 0.5)
    // 0 means that the metric is ignored and just the other used (fastest or smallest)
    double mode_val = 1; // 1 is default mode avg, harmonic mean between ratios
    if ( mode == "smaller" || mode == "faster" )
        mode_val = 0.5;
    else if ( mode == "smallest" || mode == "fastest" )
        mode_val = 0;

    // used either for fp or filter_size_bits, depending on what is provided
    double var_val = 1;
    // used for bin ratios
    double bins_val = 1;
    if ( mode == "smaller" || mode == "smallest" )
        var_val = mode_val;
    else if ( mode == "faster" || mode == "fastest" )
        bins_val = mode_val;

    double min_avg = 0;
    for ( auto const& params : simulations )
    {
        double var_ratio = 0;
        if ( filter_size )
            var_ratio = params.fp / min_fp;
        else
            var_ratio = params.filter_size_bits / static_cast< double >( min_filter_size );

        double const bins_ratio = params.n_bins / static_cast< double >( min_bins );
        double const avg        = ( 1 + std::pow( mode_val, 2 ) )
                           * ( ( var_ratio * bins_ratio ) / ( ( var_val * var_ratio ) + ( bins_val * bins_ratio ) ) );


        if ( avg < min_avg || min_avg == 0 )
        {
            min_avg = avg;
            if ( filter_size )
            {
                ibf_config.bin_size_bits =
                    ( filter_size / static_cast< double >( optimal_bins( params.n_bins ) ) ) * 8388608u;
                ibf_config.max_fp = params.fp;
            }
            else
            {
                ibf_config.bin_size_bits = params.filter_size_bits / optimal_bins( params.n_bins );
                ibf_config.max_fp        = max_fp;
            }

            ibf_config.max_hashes_bin = params.n_hashes;
            ibf_config.n_bins         = params.n_bins;
            ibf_config.hash_functions = get_optimal_hash_functions(
                ibf_config.bin_size_bits, params.n_hashes, hash_functions, max_hash_functions );
        }
    }
}


TBinMapHash create_bin_map_hash( IBFConfig const& ibf_config, THashesCount const& hashes_count )
{

    /*
     * creates binnos based on the IBF parameters and current hashes
     * distribute hashes by average on split bins
     * outputs a map with binno and indices of the hashes array to be inserted to the IBF
     */
    uint64_t    binno = 0;
    TBinMapHash bin_map_hash;
    for ( auto const& [target, count] : hashes_count )
    {
        // average hashes target
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( ibf_config.max_hashes_bin ) );
        uint64_t n_hashes_bin  = std::ceil( count / static_cast< double >( n_bins_target ) );

        if ( n_hashes_bin > ibf_config.max_hashes_bin )
            n_hashes_bin = ibf_config.max_hashes_bin;

        for ( uint64_t i = 0; i < n_bins_target; ++i )
        {
            uint64_t hashes_idx_st = i * n_hashes_bin;
            uint64_t hashes_idx_en = hashes_idx_st + n_hashes_bin - 1;
            if ( hashes_idx_st >= count )
                break;
            if ( hashes_idx_en >= count )
                hashes_idx_en = count - 1;
            bin_map_hash[binno] = std::make_tuple( target, hashes_idx_st, hashes_idx_en );
            binno++;
        }
    }
    // assert if number of bins created are same as expected
    assert( bin_map_hash.size() == ibf_config.n_bins );
    return bin_map_hash;
}

void build( TIBF&                       ibf,
            std::atomic< std::size_t >& bin_batches,
            const uint64_t              max_batch,
            const uint64_t              batch_size,
            const TBinMapHash&          bin_map_hash,
            std::string                 tmp_output_folder )
{
    /*
     * build the IBF from minimizers files and a previously generated map
     */
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

        // store files and the hashes into a map for quick access
        // since the bin batch could have repeated and muliple files
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

void print_stats( Stats& stats, const IBFConfig& ibf_config, const StopClock& timeGanonBuild )
{
    /*
     * print general statistic of the build
     */
    double elapsed = timeGanonBuild.elapsed();
    std::cerr << "ganon-build processed " << stats.total.sequences << " sequences / " << stats.total.files << " files ("
              << stats.total.length_bp / 1000000.0 << " Mbp) in " << elapsed << " seconds ("
              << ( stats.total.length_bp / 1000000.0 ) / ( elapsed / 60.0 ) << " Mbp/m)" << std::endl;

    if ( stats.total.invalid_files > 0 )
        std::cerr << " - " << stats.total.invalid_files << " invalid files skipped" << std::endl;
    if ( stats.total.skipped_sequences > 0 )
        std::cerr << " - " << stats.total.skipped_sequences << " sequences skipped" << std::endl;

    std::cerr << std::fixed << std::setprecision( 4 ) << " - max. false positive: " << ibf_config.true_max_fp;
    std::cerr << std::fixed << std::setprecision( 4 ) << " (avg.: " << ibf_config.true_avg_fp << ")" << std::endl;
    std::cerr << std::fixed << std::setprecision( 2 ) << " - filter size: "
              << ( optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits ) / static_cast< double >( 8388608u )
              << "MB" << std::endl;
}

void print_stats_verbose( const StopClock& timeGanonBuild,
                          const StopClock& timeCountStoreHashes,
                          const StopClock& timeEstimateParams,
                          const StopClock& timeBuildIBF,
                          const StopClock& timeWriteIBF )
{
    /*
     * print advanced statistic and times of the build
     */
    using ::operator<<;
    std::cerr << "Count/save hashes start: " << StopClock_datetime( timeCountStoreHashes.begin() ) << std::endl;
    std::cerr << "                    end: " << StopClock_datetime( timeCountStoreHashes.end() ) << std::endl;
    std::cerr << "            elapsed (s): " << timeCountStoreHashes.elapsed() << std::endl;
    std::cerr << "Estimate params   start: " << StopClock_datetime( timeEstimateParams.begin() ) << std::endl;
    std::cerr << "                    end: " << StopClock_datetime( timeEstimateParams.end() ) << std::endl;
    std::cerr << "            elapsed (s): " << timeEstimateParams.elapsed() << std::endl;
    std::cerr << "Building filter   start: " << StopClock_datetime( timeBuildIBF.begin() ) << std::endl;
    std::cerr << "                    end: " << StopClock_datetime( timeBuildIBF.end() ) << std::endl;
    std::cerr << "            elapsed (s): " << timeBuildIBF.elapsed() << std::endl;
    std::cerr << "Saving filer      start: " << StopClock_datetime( timeWriteIBF.begin() ) << std::endl;
    std::cerr << "                    end: " << StopClock_datetime( timeWriteIBF.end() ) << std::endl;
    std::cerr << "            elapsed (s): " << timeWriteIBF.elapsed() << std::endl;
    std::cerr << "ganon-build       start: " << StopClock_datetime( timeGanonBuild.begin() ) << std::endl;
    std::cerr << "                    end: " << StopClock_datetime( timeGanonBuild.end() ) << std::endl;
    std::cerr << "            elapsed (s): " << timeGanonBuild.elapsed() << std::endl;
    std::cerr << std::endl;
}

} // namespace detail

bool run( Config config )
{

    // Validate user provided arguments
    if ( !config.validate() )
        return false;

    // Print arguments if verbose is active
    if ( config.verbose )
        std::cerr << config;

    // Start timer for total build time
    StopClock timeGanonBuild;
    StopClock timeCountStoreHashes;
    StopClock timeEstimateParams;
    StopClock timeBuildIBF;
    StopClock timeWriteIBF;
    timeGanonBuild.start();

    // Init. Stats (general) and Totals (one for each thread)
    detail::Stats                stats;
    std::vector< detail::Total > totals( config.threads );

    // create IBF configuration and set-up fixed parameters
    IBFConfig ibf_config;
    ibf_config.kmer_size   = config.kmer_size;
    ibf_config.window_size = config.window_size;

    // Map to store number of hashes for each target {target: count}
    detail::THashesCount hashes_count;

    // Parse valid input file into a queue of InputFileMap by target and initialize hashes_count for each target or seqid
    SafeQueue< detail::InputFileMap > ifm_queue;
    for ( auto const& [target, files] :
          detail::parse_input_file( config.input_file, hashes_count, config.quiet, stats ) )
    {
        // Add to Target and Files to the SafeQueue
        ifm_queue.push( detail::InputFileMap{ target, files } );
    }
    ifm_queue.notify_push_over();

    if ( ifm_queue.empty() )
    {
        std::cerr << "No valid input files" << std::endl;
        return false;
    }

    // Create temporary output folder if not existing to write minimizer hashes
    if ( config.tmp_output_folder != "" && !std::filesystem::exists( config.tmp_output_folder ) )
    {
        std::filesystem::create_directory( config.tmp_output_folder );
    }
    else
    {
        // Delete .min hashes files in case they were previously created
        detail::delete_hashes( hashes_count, config.tmp_output_folder );
    }

    // Initialize in parallel (by target) the hash counting and storing
    timeCountStoreHashes.start();
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
    // Wait until all threads are over
    for ( auto&& task : tasks_count )
    {
        task.get();
    }
    timeCountStoreHashes.stop();

    // Sum totals from threads
    stats.add_totals( totals );

    // Define optimal parameters based on provided filter size or max.fp
    // fills ibf_config with optimal value
    timeEstimateParams.start();
    // Optimal max hashes per bin based on filter size (smallest harm.mean between fp and n. bins)
    detail::optimal_hashes( config.max_fp,
                            config.filter_size,
                            ibf_config,
                            hashes_count,
                            config.hash_functions,
                            config.max_hash_functions,
                            config.mode );
    // Calculate true final fp and average based on selected ibf_config params
    std::tie( ibf_config.true_max_fp, ibf_config.true_avg_fp ) = detail::true_false_positive(
        hashes_count, ibf_config.max_hashes_bin, ibf_config.bin_size_bits, ibf_config.hash_functions );
    timeEstimateParams.stop();

    // Print verbose arguments for ibf
    if ( config.verbose )
    {
        std::cerr << ibf_config;
        std::cerr << "Filter size: " << ( detail::optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits )
                  << " Bits";
        std::cerr << " ("
                  << ( detail::optimal_bins( ibf_config.n_bins ) * ibf_config.bin_size_bits )
                         / static_cast< double >( 8388608u )
                  << " Megabytes)" << std::endl;
    }

    if ( ibf_config.n_bins == 0 )
    {
        std::cerr << "No valid sequences to build" << std::endl;
        return false;
    }

    // Split hashes into optimal size creating technical bins
    // {binno: (target, idx_hashes_start, idx_hashes_end)}
    const detail::TBinMapHash bin_map_hash = detail::create_bin_map_hash( ibf_config, hashes_count );

    // Create IBF
    timeBuildIBF.start();
    auto ibf = detail::TIBF{ seqan3::bin_count{ ibf_config.n_bins },
                             seqan3::bin_size{ ibf_config.bin_size_bits },
                             seqan3::hash_function_count{ ibf_config.hash_functions } };

    // Add minimizer hashes to the IBF reading from previous written .min files
    // parallelize for every batch_size (64) to guarantee thread safety to the IBF
    uint64_t batch_size = 64;
    // calculate max number of batches necessary
    uint64_t max_batch = std::ceil( bin_map_hash.size() / static_cast< double >( batch_size ) );
    // Keep track of batches processed in the threads with an atomic integer
    std::atomic< std::size_t >         bin_batches = 0;
    std::vector< std::future< void > > tasks_build;
    for ( uint16_t taskn = 0; taskn < config.threads; ++taskn )
    {
        tasks_build.emplace_back( std::async( std::launch::async, [=, &ibf, &bin_batches, &bin_map_hash]() {
            detail::build( ibf, bin_batches, max_batch, batch_size, bin_map_hash, config.tmp_output_folder );
        } ) );
    }
    // Wait until all threads are over
    for ( auto&& task : tasks_build )
    {
        task.get();
    }
    timeBuildIBF.stop();

    // Delete .min hashes files
    detail::delete_hashes( hashes_count, config.tmp_output_folder );

    // Write IBF and other data structures to file
    timeWriteIBF.start();
    detail::save_filter( config, ibf, ibf_config, hashes_count, bin_map_hash );
    timeWriteIBF.stop();

    // Stop timer for total build time
    timeGanonBuild.stop();

    // Print stats and times
    if ( !config.quiet )
    {
        if ( config.verbose )
        {
            detail::print_stats_verbose(
                timeGanonBuild, timeCountStoreHashes, timeEstimateParams, timeBuildIBF, timeWriteIBF );
        }
        detail::print_stats( stats, ibf_config, timeGanonBuild );
    }

    return true;
}

} // namespace GanonBuild
