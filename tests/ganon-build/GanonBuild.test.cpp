#include "aux/Aux.hpp"

#include <cereal/archives/binary.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/std/ranges>

#include <iostream>

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>
#include <utils/adjust_seed.hpp>

#include <catch2/catch.hpp>

using namespace seqan3::literals;


namespace config_build
{

void validate_filter( const GanonBuild::Config cfg )
{
    // validate properties of the filter

    // check if file exists
    REQUIRE( std::filesystem::exists( cfg.output_file ) );
    // file should not be empty
    REQUIRE_FALSE( aux::fileIsEmpty( cfg.output_file ) );
    // load filter
    IBFConfig                                          ibf_config;
    std::vector< std::tuple< uint64_t, std::string > > bin_map;
    auto                                               filter = aux::load_ibf( cfg.output_file, ibf_config, bin_map );

    // check bin count
    REQUIRE( filter.bin_count() == bin_map.size() );

    // check hash functions
    REQUIRE( filter.hash_function_count() == ibf_config.hash_functions );
}
/*
void validate_elements( const GanonBuild::Config cfg, const sequences_type& seqs)
{
    // check if elements were properly inserted in the IBF
    // expects unique minimizers among all sequences without errors
    auto filter = aux::load_ibf( cfg.output_file );
    auto agent  = filter.counting_agent< uint16_t >();
    auto minimiser_hash = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ cfg.kmer_size } },
                                                            seqan3::window_size{ cfg.window_size },
                                                            seqan3::seed{ raptor::adjust_seed( cfg.kmer_size ) } );

    std::vector< uint16_t >             expected_counts( filter.bin_count(), 0 );
    seqan3::counting_vector< uint16_t > counts( filter.bin_count(), 0 );

    int i = 0;
    for ( auto& seq : seqs )
    {
        std::vector< uint64_t > hashes = seq | minimiser_hash | seqan3::views::to< std::vector >;
        counts += agent.bulk_count( hashes );
        // Calculate expected number of subsequences to be found (no errors==all hashes)
        expected_counts[bins[i]] += hashes.size();
        i += 1;
    }

    for ( size_t i = 0; i < counts.size(); ++i )
    {
        REQUIRE(counts[i] < expected_counts[i])
    }

}
*/


GanonBuild::Config defaultConfig( const std::string prefix )
{
    GanonBuild::Config cfg;
    cfg.output_file    = prefix + ".ibf";
    cfg.verbose        = true;
    cfg.quiet          = false;
    cfg.kmer_size      = 19;
    cfg.window_size    = 32;
    cfg.hash_functions = 0;
    cfg.max_fp         = 0.05;
    return cfg;
}


GanonBuild::Config defaultConfig( const std::string prefix, const sequences_type& seqs )
{

    // Make config
    GanonBuild::Config cfg   = defaultConfig( prefix );
    auto               files = aux::write_sequences_files( prefix, "fasta", seqs );
    aux::write_input_file( prefix + "_input.tsv", files );
    cfg.input_file = prefix + "_input.tsv";
    return cfg;
}

} // namespace config_build


// Default sequences to build
const sequences_type seqs{ "TTCAATTCGGCGTACTCAGCATCGCAGCTAGCTGTACGGCTAGTCGTCAT"_dna4,
                           "TTGGGGCTAAACAGCACTATACAGGCGGCTAGCATGTATTAGGGGAGCTC"_dna4,
                           "ACCTTCGATTTCTTTAGATCGGGGATGATGATGCATGATGCTTAGGGATT"_dna4 };

SCENARIO( "building indices", "[ganon-build]" )
{

    std::string folder_prefix{ "ganon-build-build/" };
    std::filesystem::create_directory( folder_prefix );

    SECTION( "--input file" )
    {

        SECTION( "one column - filename as target" )
        {
            auto cfg = config_build::defaultConfig( folder_prefix + "input_file_one_col", seqs );
            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            // config_build::validate_elements( cfg, seqs );
        }
    }
}
