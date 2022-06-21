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
#include <robin_hood.h>


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

    // check if bin count from filter matches the map and configuration
    REQUIRE( filter.bin_count() == bin_map.size() );
    REQUIRE( filter.bin_count() == ibf_config.n_bins );

    // check if hash functions from filter matches the configuration
    REQUIRE( filter.hash_function_count() == ibf_config.hash_functions );
}


void validate_elements( const GanonBuild::Config cfg, const SeqTarget& seqtarget )
{
    // check if elements were properly inserted in the IBF
    // expects unique minimizers among all sequences without errors
    IBFConfig                                          ibf_config;
    std::vector< std::tuple< uint64_t, std::string > > bin_map;
    auto                                               filter = aux::load_ibf( cfg.output_file, ibf_config, bin_map );
    // create map with targets and binnos
    robin_hood::unordered_map< std::string, std::vector< uint64_t > > targets;
    for ( auto const& [binno, target] : bin_map )
    {
        targets[target].push_back( binno );
    }

    seqan3::debug_stream << seqtarget.targets << std::endl;
    auto agent          = filter.counting_agent< uint16_t >();
    auto minimiser_hash = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ cfg.kmer_size } },
                                                         seqan3::window_size{ cfg.window_size },
                                                         seqan3::seed{ raptor::adjust_seed( cfg.kmer_size ) } );

    int i = 0;
    for ( auto& seq : seqtarget.sequences )
    {
        auto hashes     = seq | minimiser_hash | seqan3::views::to< std::vector >;
        auto counts_seq = agent.bulk_count( hashes );

        seqan3::debug_stream << counts_seq << std::endl;
        uint16_t expected_counts = 0;

        // sum matches for all bins of the targets this sequences belongs
        for ( auto binno : targets[seqtarget.targets[i]] )
        {
            expected_counts += counts_seq[binno];
        }

        REQUIRE( expected_counts == hashes.size() );
        i += 1;
    }
}


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


GanonBuild::Config defaultConfig( const std::string prefix, const SeqTarget& seqtarget, bool write_target )
{

    // Make config
    GanonBuild::Config cfg = defaultConfig( prefix );
    aux::write_sequences_files( seqtarget );
    cfg.input_file = aux::write_input_file( prefix, seqtarget, write_target );
    return cfg;
}

} // namespace config_build


// Default sequences to build
sequences_type seqs{ "TTCAATTCGGCGTACTCAGCATCGCAGCTAGCTGTACGGCTAGTCGTCAT"_dna4,
                     "TTGGGGCTAAACAGCACTATACAGGCGGCTAGCATGTATTAGGGGAGCTC"_dna4,
                     "ACCTTCGATTTCTTTAGATCGGGGATGATGATGCATGATGCTTAGGGATT"_dna4 };

SCENARIO( "building indices", "[ganon-build]" )
{

    std::string folder_prefix{ "ganon-build-build/" };
    std::filesystem::create_directory( folder_prefix );

    SECTION( "--input-file" )
    {

        SECTION( "one column - filename as target" )
        {
            std::string prefix{ folder_prefix + "input_file_one_col" };
            auto        seqtarget = SeqTarget( prefix, seqs );

            auto cfg = config_build::defaultConfig( prefix, seqtarget, false );

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "two columns - custom target" )
        {
            std::string                prefix{ folder_prefix + "input_file_two_cols" };
            std::vector< std::string > targets{ "T1", "T1", "T2" };
            auto                       seqtarget = SeqTarget( prefix, seqs, targets );
            auto                       cfg       = config_build::defaultConfig( prefix, seqtarget, true );

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }
    }
}
