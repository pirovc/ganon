#include "aux/Aux.hpp"

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/std/ranges>
#include <utils/dna4_traits.hpp>

#include <iostream>

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>
#include <utils/adjust_seed.hpp>

#include <catch2/catch.hpp>
#include <robin_hood.h>


using namespace seqan3::literals;


namespace config_build
{

using sequences_type       = std::vector< seqan3::dna4_vector >;
using sequence_record_type = seqan3::sequence_record< seqan3::type_list< std::vector< seqan3::dna4 >, std::string >,
                                                      seqan3::fields< seqan3::field::seq, seqan3::field::id > >;

struct SeqTarget
{
    SeqTarget()
    {
    }

    SeqTarget( std::string prefix, sequences_type& _sequences )
    {
        sequences = _sequences;
        for ( size_t i = 0; i < sequences.size(); i++ )
        {
            auto seqid  = "SEQ" + std::to_string( i );
            auto p      = std::filesystem::path( prefix ).parent_path().string();
            auto f      = std::filesystem::path( prefix ).filename().string();
            auto suffix = "." + seqid + ".fasta";
            files.push_back( std::filesystem::canonical( p ).string() + "/" + f + suffix );
            targets.push_back( f + suffix );
            headers.push_back( seqid );
        }
    }

    SeqTarget( std::string prefix, sequences_type& _sequences, std::vector< std::string >& _targets )
    {
        sequences      = _sequences;
        targets        = _targets;
        custom_targets = true;
        for ( size_t i = 0; i < sequences.size(); i++ )
        {
            auto seqid  = "SEQ" + std::to_string( i );
            auto p      = std::filesystem::path( prefix ).parent_path().string();
            auto f      = std::filesystem::path( prefix ).filename().string();
            auto suffix = "." + seqid + ".fasta";
            files.push_back( std::filesystem::canonical( p ).string() + "/" + f + suffix );
            headers.push_back( seqid );
        }
    }

    void write_input_file( const std::string input_file )
    {
        std::ofstream output_file{ input_file };
        for ( uint16_t i = 0; i < files.size(); ++i )
        {
            output_file << files[i];
            if ( custom_targets || sequence_as_target )
            {
                output_file << '\t' << targets[i];
                if ( sequence_as_target )
                {
                    output_file << '\t' << headers[i];
                }
            }

            output_file << '\n';
        }
        output_file.close();
    }

    void write_sequences_files()
    {
        for ( uint16_t i = 0; i < files.size(); ++i )
        {
            seqan3::sequence_file_output fout{ files[i] };
            sequence_record_type         rec{ sequences[i], headers[i] };
            fout.push_back( rec );
        }
    }

    std::vector< std::string >         files;
    std::vector< std::string >         headers;
    std::vector< seqan3::dna4_vector > sequences;
    std::vector< std::string >         targets;
    bool                               sequence_as_target = false;
    bool                               custom_targets     = false;
};


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
    if ( cfg.hash_functions > 0 )
        REQUIRE( filter.hash_function_count() == cfg.hash_functions );

    if ( !cfg.filter_size )
    {
        // verify if requested max fp and avg was achieved (or lower) to two decimal places
        REQUIRE( ( std::floor( ibf_config.true_max_fp * 100.0 ) / 100.0 )
                 <= ( std::floor( cfg.max_fp * 100.0 ) / 100.0 ) );
        REQUIRE( ( std::floor( ibf_config.true_avg_fp * 100.0 ) / 100.0 )
                 <= ( std::floor( cfg.max_fp * 100.0 ) / 100.0 ) );
    }
}


void validate_elements( const GanonBuild::Config cfg, const SeqTarget& seqtarget )
{
    /*
     * check if elements were properly inserted in the IBF by querying them against the generated IBF
     */

    // load IBF and data
    IBFConfig                                          ibf_config;
    std::vector< std::tuple< uint64_t, std::string > > bin_map;
    auto                                               filter = aux::load_ibf( cfg.output_file, ibf_config, bin_map );

    // agent to query and hash
    auto agent          = filter.counting_agent< uint16_t >();
    auto minimiser_hash = seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ cfg.kmer_size } },
                                                         seqan3::window_size{ cfg.window_size },
                                                         seqan3::seed{ raptor::adjust_seed( cfg.kmer_size ) } );

    // create map with targets and respective binnos {target: (binnos)}
    robin_hood::unordered_map< std::string, std::vector< uint64_t > > targets;
    for ( auto const& [binno, target] : bin_map )
    {
        targets[target].push_back( binno );
    }

    int i = 0;
    // For each sequence used to build the filter
    for ( auto& seq : seqtarget.sequences )
    {
        auto hashes     = seq | minimiser_hash | seqan3::views::to< std::vector >;
        auto counts_seq = agent.bulk_count( hashes );

        // get bins from target of this sequence
        std::vector< uint64_t > target_bins;
        if ( seqtarget.sequence_as_target )
            target_bins = targets[seqtarget.headers[i]];
        else
            target_bins = targets[seqtarget.targets[i]];

        // sum matches for all bins of the targets this sequences belongs
        uint16_t expected_counts = 0;
        for ( auto binno : target_bins )
        {
            expected_counts += counts_seq[binno];
        }

        // it should have all hashes matching among the respective bins
        REQUIRE( expected_counts == hashes.size() );
        i += 1;
    }
}


GanonBuild::Config defaultConfig( const std::string prefix )
{
    GanonBuild::Config cfg;
    cfg.input_file     = prefix + "_input.tsv";
    cfg.output_file    = prefix + ".ibf";
    cfg.verbose        = false;
    cfg.quiet          = true;
    cfg.kmer_size      = 19;
    cfg.window_size    = 32;
    cfg.hash_functions = 0;
    cfg.max_fp         = 0.05;
    return cfg;
}

} // namespace config_build

// Default sequences to build

config_build::sequences_type seqs{
    "ACACTCTTTGAAAATGCATATAATATTGAACGTTATTTTGAAATAGATTAATTACTCATATCCATTTGCTAATCTTATCG"_dna4,
    "TTTATTATATGTAATTATAAATTTATCGTTAAGCTTGACATAAGTGAGTGTATCTATGTTCTTAACAAATACATCGCGTT"_dna4,
    "TTTTATTTTTATTTCTTATGCACAAGAATAAATTATATGCATATGATAATTTCTCATTCAATGCGGATGTACATTATGGT"_dna4,
    "TATGGTAAGCTATTATGGCATGATAAAAAACCAGTCATATACCCATTGGCATCCTTATCTGATTATACTTATTATAACGA"_dna4,
    "ATCCGACCCATTTGAAACGATTTATTATGTGGAGCAATACTATAAAATTAGCTTAAATGAGAGTAAGCGAATTCAAGAAC"_dna4,
    "AAAAGGACATTTACGCACACCTTCAATTAAAACATAATAAATCATTAATTACAGCAAATGTAACGTTACATAATAAAAGT"_dna4,
    "AATAGTTCGTATTATGTTCATCGGATGAATTTACCAGCAAACATCCATGAATCACCTTACTCTCCTTTGTGCAGTGGTTC"_dna4,
    "TTTTTTAATCGTAACAAATAACATACGGTTAGATTATATAAGAAAAATTACATGCCGATTTGATTTGTGGATAAAAAAAT"_dna4,
    "CTGACTGGATAGAAATATCACCCGGAGAAAAACTCTCATACACAGTAAATTTGAATGACTATTATGCTTTTCTCCCTGCG"_dna4,
    "ATGCATCAATATGATATAGGAACTGTAGAGTTCACATTGGTAAATAGTAATTGGTTCTTAGAACAGCATATTTATGATCT"_dna4
};

SCENARIO( "building indices", "[ganon-build]" )
{

    std::string folder_prefix{ "ganon-build/" };
    std::filesystem::create_directory( folder_prefix );

    SECTION( "--input-file" )
    {

        SECTION( "one column - filename as default target" )
        {
            std::string prefix{ folder_prefix + "input_file_one_col" };
            auto        cfg = config_build::defaultConfig( prefix );

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "two columns - custom target" )
        {
            std::string prefix{ folder_prefix + "input_file_two_cols" };
            auto        cfg = config_build::defaultConfig( prefix );

            std::vector< std::string > targets{ "T1", "T9", "T1", "T8", "T1", "T1", "T1", "T1", "T4", "T1" };

            auto seqtarget = config_build::SeqTarget( prefix, seqs, targets );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "three columns - sequence as target" )
        {
            std::string prefix{ folder_prefix + "input_file_three_cols" };
            auto        cfg = config_build::defaultConfig( prefix );

            auto seqtarget               = config_build::SeqTarget( prefix, seqs );
            seqtarget.sequence_as_target = true;
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }
    }

    SECTION( "--max-fp" )
    {
        std::string prefix{ folder_prefix + "max_fp_0.01" };
        auto        cfg = config_build::defaultConfig( prefix );
        cfg.max_fp      = 0.01;

        auto seqtarget = config_build::SeqTarget( prefix, seqs );
        seqtarget.write_input_file( cfg.input_file );
        seqtarget.write_sequences_files();

        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg );
        config_build::validate_elements( cfg, seqtarget );


        std::string prefix2{ folder_prefix + "max_fp_0.5" };
        auto        cfg2 = config_build::defaultConfig( prefix2 );
        cfg2.max_fp      = 0.5;

        auto seqtarget2 = config_build::SeqTarget( prefix2, seqs );
        seqtarget2.write_input_file( cfg2.input_file );
        seqtarget2.write_sequences_files();

        REQUIRE( GanonBuild::run( cfg2 ) );
        config_build::validate_filter( cfg2 );
        config_build::validate_elements( cfg2, seqtarget2 );

        // higher max_fp should generate smaller filter
        REQUIRE( aux::fileSizeBytes( cfg.output_file ) > aux::fileSizeBytes( cfg2.output_file ) );
    }

    SECTION( "--filter-size" )
    {

        std::string prefix{ folder_prefix + "filter_size_0.1" };
        auto        cfg = config_build::defaultConfig( prefix );
        cfg.filter_size = 0.1;

        auto seqtarget = config_build::SeqTarget( prefix, seqs );
        seqtarget.write_input_file( cfg.input_file );
        seqtarget.write_sequences_files();

        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg );
        config_build::validate_elements( cfg, seqtarget );


        std::string prefix2{ folder_prefix + "filter_size_1" };
        auto        cfg2 = config_build::defaultConfig( prefix2 );
        cfg2.filter_size = 1;

        auto seqtarget2 = config_build::SeqTarget( prefix2, seqs );
        seqtarget2.write_input_file( cfg2.input_file );
        seqtarget2.write_sequences_files();

        REQUIRE( GanonBuild::run( cfg2 ) );
        config_build::validate_filter( cfg2 );
        config_build::validate_elements( cfg2, seqtarget2 );

        // check if changed filter size
        REQUIRE( aux::fileSizeBytes( cfg.output_file ) < aux::fileSizeBytes( cfg2.output_file ) );
    }

    SECTION( "--hash-functions" )
    {

        SECTION( "0 auto-detect" )
        {
            std::string prefix{ folder_prefix + "hash_functions_0" };
            auto        cfg    = config_build::defaultConfig( prefix );
            cfg.hash_functions = 0;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "2 fixed" )
        {
            std::string prefix{ folder_prefix + "hash_functions_2" };
            auto        cfg    = config_build::defaultConfig( prefix );
            cfg.hash_functions = 2;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }


        SECTION( "6 invalid" )
        {
            std::string prefix{ folder_prefix + "hash_functions_6" };
            auto        cfg    = config_build::defaultConfig( prefix );
            cfg.hash_functions = 6;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE_FALSE( GanonBuild::run( cfg ) );
        }
    }

    SECTION( "--window-size, --kmer-size" )
    {

        SECTION( "w32, k19" )
        {
            std::string prefix{ folder_prefix + "hash_functions_w32_k19" };
            auto        cfg = config_build::defaultConfig( prefix );
            cfg.window_size = 32;
            cfg.kmer_size   = 19;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "w23, k21" )
        {
            std::string prefix{ folder_prefix + "hash_functions_w23_k21" };
            auto        cfg = config_build::defaultConfig( prefix );
            cfg.window_size = 23;
            cfg.kmer_size   = 21;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "w27, k27" )
        {
            std::string prefix{ folder_prefix + "hash_functions_w27_k27" };
            auto        cfg = config_build::defaultConfig( prefix );
            cfg.window_size = 27;
            cfg.kmer_size   = 27;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "w42, k35 invalid - kmer too large" )
        {
            std::string prefix{ folder_prefix + "hash_functions_w42_k35" };
            auto        cfg = config_build::defaultConfig( prefix );
            cfg.window_size = 42;
            cfg.kmer_size   = 35;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE_FALSE( GanonBuild::run( cfg ) );
        }

        SECTION( "w12, k32 invalid - window larger than kmer" )
        {
            std::string prefix{ folder_prefix + "hash_functions_w12_k32" };
            auto        cfg = config_build::defaultConfig( prefix );
            cfg.window_size = 12;
            cfg.kmer_size   = 32;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE_FALSE( GanonBuild::run( cfg ) );
        }
    }

    SECTION( "--tmp-output-folder" )
    {

        SECTION( "empty" )
        {
            std::string prefix{ folder_prefix + "tmp_output_folder_empty" };
            auto        cfg       = config_build::defaultConfig( prefix );
            cfg.tmp_output_folder = "";

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "existing" )
        {

            std::string prefix{ folder_prefix + "tmp_output_folder_existing" };
            auto        cfg       = config_build::defaultConfig( prefix );
            cfg.tmp_output_folder = prefix;

            std::filesystem::create_directory( prefix );

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE( GanonBuild::run( cfg ) );
            config_build::validate_filter( cfg );
            config_build::validate_elements( cfg, seqtarget );
        }

        SECTION( "non-existing" )
        {

            std::string prefix{ folder_prefix + "tmp_output_folder_non_existing" };
            auto        cfg       = config_build::defaultConfig( prefix );
            cfg.tmp_output_folder = prefix;

            auto seqtarget = config_build::SeqTarget( prefix, seqs );
            seqtarget.write_input_file( cfg.input_file );
            seqtarget.write_sequences_files();

            REQUIRE_FALSE( GanonBuild::run( cfg ) );
        }
    }
}
