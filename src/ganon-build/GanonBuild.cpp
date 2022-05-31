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
    return -static_cast< uint8_t >(std::log(max_fp));
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

inline uint64_t optimal_bins( uint64_t nbins )
{
    return std::ceil( nbins / 64.0 ) * 64;
}

inline double false_positive( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t elements )
{
    return std::pow( 1 - std::exp( -hash_functions / ( bin_size_bits / static_cast< double >( elements ) ) ),
                     hash_functions );
}

std::tuple< uint64_t, uint8_t > optimal_elements_size( THashesCount const& hashes_count,
                               uint64_t            min_hashes,
                               uint64_t            max_hashes,
                               uint64_t            filter_size,
                               uint8_t             hash_functions,
                               uint8_t             max_hash_functions )
{

    
    double min_fp = 1;
    uint64_t optimal_elements = max_hashes;
    
    int optimal_hash_functions = hash_functions;
    // total + 1 zero index
    for (size_t n = max_hashes + 1; n > min_hashes; n-=100 ){
        uint64_t n_hashes = n - 1;
        uint64_t nbins = optimal_bins(n_bins(hashes_count, n_hashes));
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );
        int64_t bin_size_bits = ( filter_size / static_cast< double >( nbins ) ) * 8388608u;
        uint8_t hf = hash_functions;
        if (hf==0)
            hf = hash_functions_from_ratio(bin_size_bits, n_hashes);
        if (hf > max_hash_functions || hf == 0)
            hf = max_hash_functions; 

        double fp = false_positive(bin_size_bits, hf, n_hashes );
        double real_fp = 1 - std::pow(1 - fp, static_cast<double>(max_split_bins));

        if (real_fp < min_fp){
            min_fp = real_fp;
            optimal_hash_functions = hf;
            optimal_elements = n_hashes;
        }
    }

    return std::make_tuple( optimal_elements, optimal_hash_functions );
}

std::tuple< uint64_t, uint8_t > optimal_elements_fp( THashesCount const& hashes_count,
                               uint64_t            min_hashes,
                               uint64_t            max_hashes,
                               double              max_fp,
                               uint8_t             hash_functions,
                               uint8_t             max_hash_functions )
{
    uint64_t min_filter_size = 0;
    uint64_t optimal_elements = max_hashes;
    
    // Define optimal hash functions based on max_fp if not provided
    uint8_t optimal_hash_functions = ( hash_functions > 0 ) ? hash_functions : hash_functions_from_fp(max_fp);
    if (optimal_hash_functions > max_hash_functions  || optimal_hash_functions == 0)
        optimal_hash_functions = max_hash_functions;

    // total + 1 zero index
    for (size_t n = max_hashes + 1; n > min_hashes; n-=100 ){

        uint64_t n_hashes = n - 1;
        uint64_t nbins = optimal_bins(n_bins(hashes_count, n_hashes));
        uint64_t max_split_bins = std::ceil( max_hashes / static_cast< double >( n_hashes ) );

        int64_t bin_size_bits = bf_size(max_fp, optimal_hash_functions, n_hashes);

        double crate = correction_rate(max_split_bins, max_fp, optimal_hash_functions);

        uint64_t filter_size_bits = bin_size_bits * nbins * crate;

        if (filter_size_bits < min_filter_size || min_filter_size == 0){
            min_filter_size = filter_size_bits;
            optimal_elements = n_hashes;
        }
    }
    return std::make_tuple( optimal_elements, optimal_hash_functions );
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

    // TODO nomenclature - n_hashes or elements

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

    uint64_t elements = max_hashes;
    uint8_t hash_functions = config.hash_functions;
    uint8_t max_hash_functions = 5;

    // Create a ibf struct
    // n_bins
    // n_hashes
    // hash_funcions
    // bin_size_bits
    // max_fp
    if (config.filter_size > 0){
        // Optimal elements based on filter size (smallest fp)
        std::tie( elements, hash_functions ) = detail::optimal_elements_size( hashes_count, min_hashes, max_hashes, config.filter_size, config.hash_functions, max_hash_functions );
    }else{
        // Optimal elements based on max_fp (smallest filter size)
        std::tie( elements, hash_functions ) = detail::optimal_elements_fp( hashes_count, min_hashes, max_hashes, config.max_fp, config.hash_functions, max_hash_functions);
    }

    std::cout << "elements per bin:" << elements << std::endl;
    std::cout << "hash_functions:" << unsigned(hash_functions) << std::endl;
    // given a number of elements per bin, create map
    uint64_t binno = 0;
    for ( auto const [t, c] : hashes_count )
    {
        uint64_t nbins = std::ceil(c/static_cast<double>(elements));
        for(uint64_t i=0; i <=nbins; ++i){
            std::cout << t << '\t' << c << '\t' << binno << std::endl; 
            binno++;
        }
       
    }
    
    // build ibf reading .min files
    // parallelize by target

    // write ibf

    // Stop timer for total build time
    timeGanon.stop();

    return true;
}

} // namespace GanonBuild
