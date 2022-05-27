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
                // target itself
                ifm[fields[0]][""] = fields[0];
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
        // no sequence as target
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
    for ( auto const [target, hash] : hashes )
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

int64_t bf_size( double max_fp, uint8_t hash_functions, uint64_t elements )
{
    return std::ceil( -static_cast< double >( elements * hash_functions )
                      / std::log( 1 - std::exp( std::log( max_fp ) / hash_functions ) ) );
}

uint8_t hash_functions_from_fp(double max_fp){
    return -static_cast< int >(std::log(max_fp));
}

uint8_t hash_functions_from_ratio(uint64_t bin_size_bits, uint64_t elements ){
    return static_cast< uint8_t >(std::log(2) * ( bin_size_bits / static_cast< double > (elements) ));
}


uint64_t n_bins(THashesCount const& hashes_count, uint64_t n_hashes){
    uint64_t nbins = 0;
    for ( auto const& [target, hcount] : hashes_count )
    {
        nbins += std::ceil( hcount / static_cast< double >( n_hashes ) );
    }
    return nbins;
}


double correction_rate( uint64_t max_split_bins, double max_fp, uint8_t hash_functions )
{
    return std::log(1 - std::exp(std::log(1.0 - std::pow(1 - max_fp, static_cast<double>(max_split_bins))) / hash_functions))
           / std::log( 1.0 - std::exp( std::log( max_fp ) / hash_functions ) );
}

inline uint64_t get_optimal_bins( uint64_t nbins )
{
    return std::ceil( nbins / 64.0 ) * 64;
}

inline double get_fp( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t elements )
{
    return std::pow( 1 - std::exp( -hash_functions / ( bin_size_bits / static_cast< double >( elements ) ) ),
                     hash_functions );
}

uint64_t optimal_elements_bin( THashesCount const& hashes_count,
                               uint64_t            min_hashes,
                               uint64_t            max_hashes,
                               uint64_t            filter_size,
                               uint8_t             hash_functions,
                               uint8_t             max_hash_functions )
{
    std::cout << "n_hashes" << '\t' << "bin_size_bits" << '\t' << "max_split_bins" << '\t' << "crate" << '\t' << "nbins" << '\t' << "fp" << std::endl;
    // total-1 zero index
    for (size_t n = max_hashes - 1; n >= min_hashes; n-=100 ){
        uint64_t n_hashes = n + 1;
        uint64_t nbins = get_optimal_bins(n_bins(hashes_count, n_hashes));
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );


        int64_t bin_size_bits = ( filter_size / static_cast< double >( nbins ) ) * 8388608u;

        if (hash_functions==0){
            hash_functions = hash_functions_from_ratio(bin_size_bits, n_hashes);
            if (hash_functions>max_hash_functions)
                hash_functions = max_hash_functions;
        }

        double fp = get_fp(bin_size_bits, hash_functions, n_hashes );
        
        
        double crate = correction_rate(max_split_bins, fp, hash_functions);

        uint64_t filter_size_bits = bin_size_bits * nbins * crate;
        std::cout << n_hashes << '\t';
        std::cout << bin_size_bits << '\t';
        std::cout << max_split_bins << '\t';
        std::cout << crate << '\t';
        std::cout << nbins << '\t';
        std::cout << fp*crate << '\t';
        std::cout << std::endl;
    }   
}

uint64_t optimal_elements_bin( THashesCount const& hashes_count,
                               uint64_t            min_hashes,
                               uint64_t            max_hashes,
                               double              max_fp,
                               uint8_t             hash_functions )
{

    std::cout << min_hashes << std::endl;
    std::cout << max_hashes << std::endl;

    std::cout << "n_hashes" << '\t' << "bin_size_bits" << '\t' << "max_split_bins" << '\t' << "crate" << '\t' << "nbins" << '\t' << "filter_size_bits" << std::endl;
        
    // total-1 zero index
    for (size_t n = max_hashes - 1; n >= min_hashes; n-=100 ){

        uint64_t n_hashes = n + 1;
        uint64_t nbins = get_optimal_bins(n_bins(hashes_count, n_hashes));
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );

        int64_t bin_size_bits = bf_size(max_fp, hash_functions, n_hashes);
        double crate = correction_rate(max_split_bins, max_fp, hash_functions);

        uint64_t filter_size_bits = bin_size_bits * nbins * crate;
        std::cout << n_hashes << '\t';
        std::cout << bin_size_bits << '\t';
        std::cout << max_split_bins << '\t';
        std::cout << crate << '\t';
        std::cout << nbins << '\t';
        std::cout << filter_size_bits << '\t';
        std::cout << std::endl;

        if (bin_size_bits < 0)
            break;

    }

    // return min filter_size
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


    // Parse input file into map {file: {"": target}} or {file: {sequence: target,...}}
    detail::TInputFilesMap ifm = detail::parse_input_file( config.input_file );

    // TODO parallelize by file
    detail::THashesCount hashes_count;
    uint64_t     min_hashes = 0;
    uint64_t     max_hashes = 0;
    for ( auto& [file, target] : ifm )
    {
        // TODO check if file exists

        detail::THashes hashes = detail::count_hashes( file, target, config );

        // store to file
        detail::store_hashes( hashes, config );

        for ( auto const [t, hash] : hashes )
        {
            hashes_count[t] = hash.size();
            //std::cout << file << '\t' << t << '\t' << hash.size() << std::endl;
            if ( hashes_count[t] > max_hashes )
                max_hashes = hashes_count[t];
            if ( hashes_count[t] < min_hashes || min_hashes == 0 )
                min_hashes = hashes_count[t];
        }
    }

    uint8_t max_hash_functions = 5;
    if (config.filter_size > 0){
        // Define hash functions if not provided
        auto elements = detail::optimal_elements_bin( hashes_count, min_hashes, max_hashes, config.filter_size, config.hash_functions, max_hash_functions );

    }else{
        // define filter size based on max_fp

        // Define hash functions if not provided
        uint8_t hash_functions = ( config.hash_functions > 0 ) ? hash_functions : detail::hash_functions_from_fp(config.max_fp);
        if (hash_functions > max_hash_functions)
            hash_functions = max_hash_functions;

        // Create technical bins -> optimize split of bins, hash functions and corrrection ratio
        auto elements = detail::optimal_elements_bin( hashes_count, min_hashes, max_hashes, config.max_fp, hash_functions );
    }

    // given a number of elements per bin, create map
    // target -> 

    // build ibf reading .min files

    // write ibf

    // Stop timer for total build time
    timeGanon.stop();

    return true;
}

} // namespace GanonBuild
