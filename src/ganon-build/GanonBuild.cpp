#include "GanonBuild.hpp"

#include <robin_hood.h>

#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>
#include <utils/adjust_seed.hpp>
#include <utils/dna4_traits.hpp>
#include <utils/load_map.hpp>

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

typedef robin_hood::unordered_map< std::string, std::string >                           TTarget;
typedef robin_hood::unordered_map< std::string, TTarget >                               TInputFilesMap;
typedef robin_hood::unordered_map< std::string, robin_hood::unordered_set< uint64_t > > THashes;
typedef robin_hood::unordered_map< std::string, uint64_t >                              THashesCount;
typedef robin_hood::unordered_map< uint64_t, std::tuple<std::string, uint64_t, uint64_t> > TBinMap;


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
};

inline std::string get_seqid( std::string header )
{
    return header.substr( 0, header.find( ' ' ) );
}

TInputFilesMap parse_input_file( const std::string& input_file )
{
    std::string    line;
    std::ifstream  infile( input_file );
    TInputFilesMap ifm;
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field ); // file [<tab> target <tab> sequence]

        if ( fields.size() == 1 )
        {
            // skip repeated files
            if ( !ifm.count( fields[0] ) )
            {
                // target itself (filename only)
                ifm[fields[0]][""] = std::filesystem::path( fields[0] ).filename();
            }
        }
        else if ( fields.size() == 2 )
        {
            // provided target
            ifm[fields[0]][""] = fields[1];
        }
        else
        {
            // sequence as target
            ifm[fields[0]][fields[2]] = fields[1];
        }
    }
    return ifm;
}

THashes count_hashes( const std::string& file, TTarget& target, const GanonBuild::Config& config )
{
    auto minimiser_view = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
                                                         seqan3::window_size{ config.window_size },
                                                         seqan3::seed{ raptor::adjust_seed( config.kmer_size ) } );

    THashes                                                                                                     hashes;
    seqan3::sequence_file_input< raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq > > fin{
        file
    };
    for ( auto& [header, seq] : fin )
    {
        const auto mh = seq | minimiser_view | std::views::common;
        // file as target
        if ( target.count( "" ) )
        {
            hashes[target[""]].insert( mh.begin(), mh.end() );
        }
        else
        {
            // hashes for each sequence
            const auto seqid = get_seqid( header );
            if ( target.count( seqid ) )
            {
                hashes[target[seqid]].insert( mh.begin(), mh.end() );
            }
        }
    }

    return hashes;
}

void store_hashes( THashes const& hashes, const GanonBuild::Config& config )
{
    // write to file
    for ( auto const& [target, hash] : hashes )
    {
        std::filesystem::path outf{ config.tmp_output_folder };
        outf += target + ".min";
        std::ofstream outfile{ outf, std::ios::binary };
        for ( auto&& h : hashes )
        {
            outfile.write( reinterpret_cast< const char* >( &h ), sizeof( h ) );
        }
    }
}

int64_t bf_size( double max_fp, uint8_t hash_functions, uint64_t n_hashes )
{
    return std::ceil( -static_cast< double >( n_hashes * hash_functions )
                      / std::log( 1 - std::exp( std::log( max_fp ) / hash_functions ) ) );
}

uint8_t hash_functions_from_fp( double max_fp )
{
    return -static_cast< uint8_t >( std::log( max_fp ) );
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
    return std::log( 1
                     - std::exp( std::log( 1.0 - std::pow( 1 - max_fp, static_cast< double >( max_split_bins ) ) )
                                 / hash_functions ) )
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

void optimal_hashes_size( IBFConfig&          ibf_config,
                          THashesCount const& hashes_count,
                          uint64_t            min_hashes,
                          uint64_t            max_hashes,
                          uint64_t            filter_size,
                          uint8_t             hash_functions,
                          uint8_t             max_hash_functions )
{


    double min_fp = 1;
    // total + 1 zero index
    for ( size_t n = max_hashes + 1; n > min_hashes; n -= 100 )
    {
        uint64_t n_hashes      = n - 1;
        uint64_t n_bins        = number_of_bins( hashes_count, n_hashes );
        int64_t  bin_size_bits = ( filter_size / static_cast< double >( optimal_bins( n_bins ) ) ) * 8388608u;
        uint8_t  hf            = hash_functions;
        if ( hf == 0 )
            hf = hash_functions_from_ratio( bin_size_bits, n_hashes );
        if ( hf > max_hash_functions || hf == 0 )
            hf = max_hash_functions;

        double   fp             = false_positive( bin_size_bits, hf, n_hashes );
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );
        double   real_fp        = 1 - std::pow( 1 - fp, static_cast< double >( max_split_bins ) );

        if ( real_fp < min_fp )
        {
            min_fp                    = real_fp;
            ibf_config.max_hashes_bin = n_hashes;
            ibf_config.bin_size_bits  = bin_size_bits;
            ibf_config.hash_functions = hf;
            ibf_config.max_fp         = real_fp;
            ibf_config.n_bins         = n_bins;
        }
    }
}

void optimal_hashes_fp( IBFConfig&          ibf_config,
                        THashesCount const& hashes_count,
                        uint64_t            min_hashes,
                        uint64_t            max_hashes,
                        double              max_fp,
                        uint8_t             hash_functions,
                        uint8_t             max_hash_functions )
{
    uint64_t min_filter_size = 0;

    // Define optimal hash functions based on max_fp if not provided
    uint8_t optimal_hash_functions = ( hash_functions > 0 ) ? hash_functions : hash_functions_from_fp( max_fp );
    if ( optimal_hash_functions > max_hash_functions || optimal_hash_functions == 0 )
        optimal_hash_functions = max_hash_functions;

    ibf_config.hash_functions = optimal_hash_functions;
    ibf_config.max_fp         = max_fp;

    // total + 1 zero index
    for ( size_t n = max_hashes + 1; n > min_hashes; n -= 100 )
    {

        uint64_t n_hashes       = n - 1;
        uint64_t n_bins         = number_of_bins( hashes_count, n_hashes );
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );

        int64_t bin_size_bits = bf_size( max_fp, optimal_hash_functions, n_hashes );

        double crate = correction_rate( max_split_bins, max_fp, optimal_hash_functions );

        uint64_t filter_size_bits = bin_size_bits * optimal_bins( n_bins ) * crate;

        if ( filter_size_bits < min_filter_size || min_filter_size == 0 )
        {
            min_filter_size           = filter_size_bits;
            ibf_config.max_hashes_bin = n_hashes;
            ibf_config.bin_size_bits  = bin_size_bits;
            ibf_config.n_bins         = n_bins;
        }
    }
}

TBinMap create_binmap(IBFConfig const&          ibf_config,
                        THashesCount const& hashes_count){

    uint64_t binno = 0;
    TBinMap binmap;
    for ( auto const& [target, count] : hashes_count )
    {
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( ibf_config.max_hashes_bin ) );
        uint64_t n_hashes_bin = std::ceil( count / static_cast< double >( n_bins_target ) );

        if ( n_hashes_bin > ibf_config.max_hashes_bin )
            n_hashes_bin = ibf_config.max_hashes_bin;

        for ( uint64_t i = 0; i < n_bins_target; ++i )
        {
            uint64_t hashes_idx_st = i * n_hashes_bin;
            uint64_t hashes_idx_en = hashes_idx_st + n_hashes_bin - 1;
            if ( hashes_idx_st >= count )
                hashes_idx_en = count - 1;
            binmap[binno] = std::make_tuple(target, hashes_idx_st, hashes_idx_en); 
            binno++;
        }
    }
    return binmap;
}

} // namespace detail

bool run( Config config )
{

    // Validate configuration input
    if ( !config.validate() )
        return false;

    if ( config.verbose )
        std::cerr << config;

    // Start timer for total build time
    StopClock timeGanon;
    timeGanon.start();

    // Limit on hash_functions from seqan3
    uint8_t max_hash_functions = 5;

    // create configuration and set up user parameters
    detail::IBFConfig ibf_config;
    ibf_config.kmer_size   = config.kmer_size;
    ibf_config.window_size = config.window_size;

    // Parse input file into map {file: {"": target or filename}} or {file: {sequence: target,...}}
    detail::TInputFilesMap ifm = detail::parse_input_file( config.input_file );

    // TODO parallelize by file
    detail::THashesCount hashes_count;
    uint64_t             min_hashes = 0;
    uint64_t             max_hashes = 0;
    for ( auto& [file, target] : ifm )
    {
        // TODO check if file exists
        detail::THashes hashes = detail::count_hashes( file, target, config );

        // store to file
        // TODO check if tmp folder exists, if not create it
        detail::store_hashes( hashes, config );

        // TODO initialize map before adding to do it in parallel
        for ( auto const [t, hash] : hashes )
        {
            hashes_count[t] = hash.size();
            // std::cout << file << '\t' << t << '\t' << hash.size() << std::endl;
            if ( hashes_count[t] > max_hashes )
                max_hashes = hashes_count[t];
            if ( hashes_count[t] < min_hashes || min_hashes == 0 )
                min_hashes = hashes_count[t];
        }
    }

    // Define optimal parameters and fill ibf_config
    if ( config.filter_size > 0 )
    {
        // Optimal max hashes per bin based on filter size (smallest fp)
        detail::optimal_hashes_size( ibf_config,
                                     hashes_count,
                                     min_hashes,
                                     max_hashes,
                                     config.filter_size,
                                     config.hash_functions,
                                     max_hash_functions );
    }
    else
    {
        // Optimal max hashes per bin based on max_fp (smallest filter size)
        detail::optimal_hashes_fp(
            ibf_config, hashes_count, min_hashes, max_hashes, config.max_fp, config.hash_functions, max_hash_functions );
    }

    std::cout << "ibf_config:" << std::endl;
    std::cout << "n_bins " << ibf_config.n_bins << std::endl;
    std::cout << "max_hashes_bin " << ibf_config.max_hashes_bin << std::endl;
    std::cout << "hash_functions " << unsigned( ibf_config.hash_functions ) << std::endl;
    std::cout << "kmer_size " << unsigned( ibf_config.kmer_size ) << std::endl;
    std::cout << "window_size " << ibf_config.window_size << std::endl;
    std::cout << "bin_size_bits " << ibf_config.bin_size_bits << std::endl;


    // Create a map {binno: (target, idx_hashes_start, idx_hashes_end)}
    detail::TBinMap binmap = create_binmap(ibf_config, hashes_count);

    // build ibf reading .min files
    // parallelize for every 64 entries in binmap, insert to ibf
    
    // write ibf and other infos

    // Stop timer for total build time
    timeGanon.stop();

    return true;
}

} // namespace GanonBuild
