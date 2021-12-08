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

#include <catch2/catch.hpp>

using namespace seqan3::literals;

// Sequences to build the ibf
const ids_type       ids{ "S1", "S2", "S3" };
const sequences_type seqs{ "TTCAATTCGGCGTACTCAGCATCGCAGCTAGCTGTACGGCTAGTCGTCAT"_dna5,
                           "TTGGGGCTAAACAGCACTATACAGGCGGCTAGCATGTATTAGGGGAGCTC"_dna5,
                           "ACCTTCGATTTCTTTAGATCGGGGATGATGATGCATGATGCTTAGGGATT"_dna5 };
const bins_type      bins{ 0, 1, 2 };

const ids_type       extra_ids{ "S4", "S5", "S6" };
const sequences_type extra_seqs{ "ATGAATTCAAGCCAATGTCGTTTGAAACAGAAGATGGTATTGCTACTGGC"_dna5,
                                 "TGCTGCCATCAACTTGCAGAAGATGTCCTTTTCTGCGGTCTACGCTCAAG"_dna5,
                                 "ATGCTGGTTAGACAGGACCTGTTAAGAAAAAGGAAACTCTCAATTGCACC"_dna5 };
const bins_type      extra_bins{ 3, 4, 5 };

namespace config_build
{

void validate_filter( const GanonBuild::Config cfg, const bins_type& bins )
{
    // validate properties of the filter

    // check if file exists
    REQUIRE( std::filesystem::exists( cfg.output_filter_file ) );
    // file should not be empty
    REQUIRE_FALSE( aux::fileIsEmpty( cfg.output_filter_file ) );
    // load filter
    seqan3::interleaved_bloom_filter<> filter = aux::load_ibf( cfg.output_filter_file );
    // check bin count
    uint32_t last_bin = *std::max_element( std::begin( bins ), std::end( bins ) );
    REQUIRE( filter.bin_count() == last_bin + 1 );
    // check hash functions
    REQUIRE( filter.hash_function_count() == cfg.hash_functions );

    // check size
    if ( cfg.filter_size_mb > 0 )
    {
        REQUIRE( filter.bit_size() == cfg.filter_size_mb * 8388608u );
    }
    else if (cfg.bin_size_bits > 0)
    {
        uint64_t optimal_bins = ( std::floor( filter.bin_count() / 64 ) + 1 ) * 64;
        REQUIRE( filter.bit_size() == cfg.bin_size_bits * optimal_bins );
    }else{
        // using --false-positive, auto calculate
        REQUIRE( filter.bit_size() > 0 );
    }
}

void validate_elements( const GanonBuild::Config cfg, const sequences_type& seqs, const bins_type& bins )
{
    // check if elements were properly inserted in the IBF
    // expects unique k-mers among all sequences without errors
    seqan3::interleaved_bloom_filter<> filter = aux::load_ibf( cfg.output_filter_file );
    auto                               agent  = filter.counting_agent();

    auto kmer_adaptor      = seqan3::views::kmer_hash( seqan3::ungapped{ cfg.kmer_size } );
    auto minimizer_adaptor = seqan3::views::minimiser_hash(
        seqan3::shape{ seqan3::ungapped{ cfg.kmer_size } }, seqan3::window_size{ cfg.window_size }, seqan3::seed{ 0 } );

    std::vector< uint16_t >             expected_output( filter.bin_count(), 0 );
    seqan3::counting_vector< uint16_t > output( filter.bin_count(), 0 );

    int i = 0;
    for ( auto& seq : seqs )
    {
        auto hashes = cfg.window_size ? seq | minimizer_adaptor | seqan3::views::to< std::vector >
                                      : seq | kmer_adaptor | seqan3::views::to< std::vector >;
        output += agent.bulk_count( hashes );
        // Calculate expected number of subsequences to be found (no errors==all hashes)
        expected_output[bins[i]] = hashes.size();
        i += 1;
    }

    // If filter was build with --false-positive rate, results may have FP counts
    if (cfg.false_positive > 0)
        for(size_t i=0; i<output.size(); ++i)
            REQUIRE( output[i] >= expected_output[i] );
    else
        REQUIRE( output == expected_output ); // no FP, find exact hashes
}


GanonBuild::Config defaultConfig( const std::string prefix )
{
    GanonBuild::Config cfg;
    cfg.bin_size_bits      = 5000;
    cfg.verbose            = false;
    cfg.quiet              = true;
    cfg.kmer_size          = 19;
    cfg.hash_functions     = 3;
    cfg.output_filter_file = prefix + ".ibf";
    return cfg;
}

GanonBuild::Config defaultConfig( const std::string     prefix,
                                  const sequences_type& seqs,
                                  const ids_type&       ids,
                                  const bins_type&      bins )
{

    // Make config
    GanonBuild::Config cfg = defaultConfig( prefix );

    // if input sequences are sent, create files
    if ( seqs.size() )
    {
        aux::write_sequences( prefix + ".fasta", seqs, ids );
        cfg.reference_files = { prefix + ".fasta" };
        aux::write_seqid_bin( prefix + "_seqid_bin.tsv", seqs, ids, bins );
        cfg.seqid_bin_file = prefix + "_seqid_bin.tsv";
    }
    return cfg;
}

} // namespace config_build


SCENARIO( "building indices", "[ganon-build]" )
{

    std::string folder_prefix{ "ganon-build-build/" };
    std::filesystem::create_directory( folder_prefix );

    SECTION( "with default conf." )
    {
        auto cfg = config_build::defaultConfig( folder_prefix + "default", seqs, ids, bins );
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with multiple --reference-files" )
    {
        std::string prefix = folder_prefix + "reference_files";
        auto        cfg    = config_build::defaultConfig( prefix );

        // write files separetly
        aux::write_sequences( prefix + ".A.fasta", seqs, ids );
        aux::write_sequences( prefix + ".B.fasta", extra_seqs, extra_ids );

        // merge entries
        auto merged_ids  = aux::vconcat( ids, extra_ids );
        auto merged_seqs = aux::vconcat( seqs, extra_seqs );
        auto merged_bins = aux::vconcat( bins, extra_bins );
        aux::write_seqid_bin( prefix + "_seqid_bin.tsv", merged_seqs, merged_ids, merged_bins );

        cfg.reference_files = { prefix + ".A.fasta", prefix + ".B.fasta" };
        cfg.seqid_bin_file  = prefix + "_seqid_bin.tsv";

        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, merged_bins );
        config_build::validate_elements( cfg, merged_seqs, merged_bins );
    }

    SECTION( "with multiple --reference-files without --seqid-bin-file" )
    {
        // without --seqid-bin-file it should create one bin per sequence file
        std::string prefix = folder_prefix + "reference_files_wo_seqid_bin";
        auto        cfg    = config_build::defaultConfig( prefix );
        // write one sequence per file
        aux::write_sequences( prefix + ".S1.fasta", { seqs[0] }, { ids[0] } );
        aux::write_sequences( prefix + ".S2.fasta", { seqs[1] }, { ids[1] } );
        aux::write_sequences( prefix + ".S3.fasta", { seqs[2] }, { ids[2] } );
        cfg.reference_files = { prefix + ".S1.fasta", prefix + ".S2.fasta", prefix + ".S3.fasta" };

        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with multiple --reference-files without --seqid-bin-file and --window-size 27 and --false-positive 0.1" )
    {
        // without --seqid-bin-file it should create one bin per sequence file
        std::string prefix = folder_prefix + "reference_files_wo_seqid_bin_window_size_27_false_positive_0.1";
        auto        cfg    = config_build::defaultConfig( prefix );
        cfg.window_size = 27;
        cfg.bin_size_bits = 0;
        cfg.false_positive = 0.1;
        // write one sequence per file
        aux::write_sequences( prefix + ".S1.fasta", { seqs[0] }, { ids[0] } );
        aux::write_sequences( prefix + ".S2.fasta", { seqs[1] }, { ids[1] } );
        aux::write_sequences( prefix + ".S3.fasta", { seqs[2] }, { ids[2] } );
        cfg.reference_files = { prefix + ".S1.fasta", prefix + ".S2.fasta", prefix + ".S3.fasta" };

        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with --window-size 23" )
    {
        auto cfg        = config_build::defaultConfig( folder_prefix + "window_size_23", seqs, ids, bins );
        cfg.window_size = 32;
        // run ganon-build
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }


    SECTION( "with --false-positive 0.05" )
    {
        auto cfg        = config_build::defaultConfig( folder_prefix + "false_positive_0.05", seqs, ids, bins );
        cfg.bin_size_bits = 0;
        cfg.false_positive = 0.5;
        // run ganon-build
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with --kmer-size 11" )
    {
        auto cfg      = config_build::defaultConfig( folder_prefix + "kmer_size_11", seqs, ids, bins );
        cfg.kmer_size = 11;
        // run ganon-build
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with --kmer-size 27" )
    {
        auto cfg      = config_build::defaultConfig( folder_prefix + "kmer_size_27", seqs, ids, bins );
        cfg.kmer_size = 27;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with --kmer-size 31" )
    {
        // Should failed due to size limitation with dna5 (max 27, otherwise 32 with dna4)
        // https://docs.seqan.de/seqan/3-master-user/group__search__views.html#ga6e598d6a021868f704d39df73252974f
        auto cfg      = config_build::defaultConfig( folder_prefix + "kmer_size_31", seqs, ids, bins );
        cfg.kmer_size = 31;
        REQUIRE_THROWS( GanonBuild::run( cfg ) );
    }

    SECTION( "with --filter-size-mb 2" )
    {
        auto cfg           = config_build::defaultConfig( folder_prefix + "filter_size_mb_2", seqs, ids, bins );
        cfg.bin_size_bits  = 0;
        cfg.filter_size_mb = 2;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with --hash-functions 2" )
    {
        auto cfg           = config_build::defaultConfig( folder_prefix + "hash_functions_2", seqs, ids, bins );
        cfg.hash_functions = 2;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }

    SECTION( "with --directory-reference-files and --extension" )
    {
        // write file with specific extension "TEST.fasta"
        auto cfg = config_build::defaultConfig( folder_prefix + "directory_reference_files.TEST", seqs, ids, bins );
        cfg.reference_files           = {};
        cfg.directory_reference_files = std::filesystem::canonical( folder_prefix );
        cfg.extension                 = ".TEST.fasta";
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, bins );
        config_build::validate_elements( cfg, seqs, bins );
    }
}


SCENARIO( "updating indices", "[ganon-build]" )
{


    std::string folder_prefix{ "ganon-build-update/" };
    std::filesystem::create_directory( folder_prefix );

    // build default filter
    auto cfg_build = config_build::defaultConfig( folder_prefix + "update_base_build", seqs, ids, bins );
    REQUIRE( GanonBuild::run( cfg_build ) );
    config_build::validate_filter( cfg_build, bins );
    config_build::validate_elements( cfg_build, seqs, bins );

    SECTION( "with --update-filter-file creating 3 new bins" )
    {
        // update it
        auto cfg_update =
            config_build::defaultConfig( folder_prefix + "3new_update", extra_seqs, extra_ids, extra_bins );
        // set filter to update
        cfg_update.update_filter_file = cfg_build.output_filter_file;
        REQUIRE( GanonBuild::run( cfg_update ) );

        auto merged_seqs = aux::vconcat( seqs, extra_seqs );
        auto merged_bins = aux::vconcat( bins, extra_bins );

        config_build::validate_filter( cfg_update, merged_bins );
        config_build::validate_elements( cfg_update, merged_seqs, merged_bins );

        SECTION( "with --update-complete" )
        {
            auto merged_ids = aux::vconcat( ids, extra_ids );
            // update complete, send all sequences
            auto cfg_update_complete = config_build::defaultConfig(
                folder_prefix + "3new_update_complete", merged_seqs, merged_ids, merged_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;
            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete, merged_bins );
            config_build::validate_elements( cfg_update_complete, merged_seqs, merged_bins );

            // should be the same as not complete
            REQUIRE( aux::filesAreEqual( cfg_update_complete.output_filter_file, cfg_update.output_filter_file ) );
        }
    }

    SECTION( "with --update-filter-file creating 3 new bins without --seqid-bin-file" )
    {
        // without --seqid-bin-file it should create one bin per sequence file
        std::string prefix            = folder_prefix + "3new_update_wo_seqid_bin";
        auto        cfg_update        = config_build::defaultConfig( prefix );
        cfg_update.update_filter_file = cfg_build.output_filter_file;

        // write one sequence per file
        aux::write_sequences( prefix + ".S4.fasta", { extra_seqs[0] }, { extra_ids[0] } );
        aux::write_sequences( prefix + ".S5.fasta", { extra_seqs[1] }, { extra_ids[1] } );
        aux::write_sequences( prefix + ".S6.fasta", { extra_seqs[2] }, { extra_ids[2] } );
        cfg_update.reference_files = { prefix + ".S4.fasta", prefix + ".S5.fasta", prefix + ".S6.fasta" };

        REQUIRE( GanonBuild::run( cfg_update ) );

        auto merged_seqs = aux::vconcat( seqs, extra_seqs );
        auto merged_bins = aux::vconcat( bins, extra_bins );

        config_build::validate_filter( cfg_update, merged_bins );
        config_build::validate_elements( cfg_update, merged_seqs, merged_bins );
    }


    SECTION( "with --update-filter-file creating 80 new bins" )
    {
        // set sequences to new bins
        const bins_type new_bins{ 12, 64, 82 };
        // update it
        auto cfg_update =
            config_build::defaultConfig( folder_prefix + "80new_update", extra_seqs, extra_ids, new_bins );
        // set filter to update
        cfg_update.update_filter_file = cfg_build.output_filter_file;
        REQUIRE( GanonBuild::run( cfg_update ) );

        auto merged_seqs = aux::vconcat( seqs, extra_seqs );
        auto merged_bins = aux::vconcat( bins, new_bins );

        config_build::validate_filter( cfg_update, merged_bins );
        config_build::validate_elements( cfg_update, merged_seqs, merged_bins );

        // new file should be bigger (double bits + overhead)
        REQUIRE( aux::fileSizeBytes( cfg_update.output_filter_file )
                 > aux::fileSizeBytes( cfg_build.output_filter_file ) );

        SECTION( "with --update-complete" )
        {
            auto merged_ids = aux::vconcat( ids, extra_ids );
            // update complete, send all sequences
            auto cfg_update_complete = config_build::defaultConfig(
                folder_prefix + "80new_update_complete", merged_seqs, merged_ids, merged_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;
            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete, merged_bins );
            config_build::validate_elements( cfg_update_complete, merged_seqs, merged_bins );

            // should be the same as not complete
            REQUIRE( aux::filesAreEqual( cfg_update_complete.output_filter_file, cfg_update.output_filter_file ) );
        }
    }

    SECTION( "with --update-filter-file without adding new bins (updating current bins)" )
    {
        // update it
        // use bins instead of extra_bins
        auto cfg_update = config_build::defaultConfig( folder_prefix + "0new_update", extra_seqs, extra_ids, bins );
        // set filter to update
        cfg_update.update_filter_file = cfg_build.output_filter_file;
        REQUIRE( GanonBuild::run( cfg_update ) );
        config_build::validate_filter( cfg_update, bins );
        // validate seqs
        config_build::validate_elements( cfg_update, seqs, bins );
        // validate extra_seqs
        config_build::validate_elements( cfg_update, extra_seqs, bins );

        SECTION( "with --update-complete" )
        {
            auto merged_seqs   = aux::vconcat( seqs, extra_seqs );
            auto merged_ids    = aux::vconcat( ids, extra_ids );
            auto repeated_bins = aux::vconcat( bins, bins );
            // update complete, send all sequences
            auto cfg_update_complete = config_build::defaultConfig(
                folder_prefix + "0new_update_complete", merged_seqs, merged_ids, repeated_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;
            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete, bins );

            // validate seqs
            config_build::validate_elements( cfg_update_complete, seqs, bins );
            // validate extra_seqs
            config_build::validate_elements( cfg_update_complete, extra_seqs, bins );

            // should be the same as not complete
            REQUIRE( aux::filesAreEqual( cfg_update_complete.output_filter_file, cfg_update.output_filter_file ) );
        }

        SECTION( "with --update-complete removing sequences" )
        {
            // create update file without first sequences

            ids_type       rem_ids{ ids[1], ids[2] };
            sequences_type rem_seqs{ seqs[1], seqs[2] };
            bins_type      rem_bins{ bins[1], bins[2] };

            auto merged_ids  = aux::vconcat( rem_ids, extra_ids );
            auto merged_seqs = aux::vconcat( rem_seqs, extra_seqs );
            auto merged_bins = aux::vconcat( rem_bins, bins );

            auto cfg_update_complete = config_build::defaultConfig(
                folder_prefix + "0new_update_complete_remove", merged_seqs, merged_ids, merged_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;

            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete, merged_bins );

            // validate seqs (without first)
            config_build::validate_elements( cfg_update_complete, rem_seqs, rem_bins );
            // validate extra_seqs (all should be there)
            config_build::validate_elements( cfg_update_complete, extra_seqs, bins );

            // should be different as before
            REQUIRE_FALSE(
                aux::filesAreEqual( cfg_update_complete.output_filter_file, cfg_update.output_filter_file ) );
        }
    }
}
