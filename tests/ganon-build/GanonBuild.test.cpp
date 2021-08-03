#include "aux/Aux.hpp"

#include <seqan3/core/debug_stream.hpp>

#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/std/ranges>

#include <filesystem>
#include <iostream>

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

using namespace seqan3::literals;

using sequences_type = std::vector< seqan3::dna5_vector >;
using ids_type       = std::vector< std::string >;
using bins_type      = std::vector< uint16_t >;

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

using sequence_record_type = seqan3::sequence_record< seqan3::type_list< std::vector< seqan3::dna5 >, std::string >,
                                                      seqan3::fields< seqan3::field::seq, seqan3::field::id > >;


void validate_filter( std::string      output_filter_file,
                      uint16_t         hash_functions,
                      uint32_t         filter_size_mb,
                      uint64_t         bin_size_bits,
                      const bins_type& bins )
{
    // validate properties of the filter

    // check if file exists
    REQUIRE( std::filesystem::exists( output_filter_file ) );
    // file should not be empty
    REQUIRE_FALSE( aux::fileIsEmpty( output_filter_file ) );
    // load filter
    seqan3::interleaved_bloom_filter<> filter = aux::load_ibf( output_filter_file );
    // check bin count
    uint32_t last_bin = *std::max_element( std::begin( bins ), std::end( bins ) );
    REQUIRE( filter.bin_count() == last_bin + 1 );
    // check hash functions
    REQUIRE( filter.hash_function_count() == hash_functions );

    // check size
    if ( filter_size_mb > 0 )
    {
        REQUIRE( filter.bit_size() == filter_size_mb * 8388608u );
    }
    else
    {
        uint64_t optimal_bins = ( std::floor( filter.bin_count() / 64 ) + 1 ) * 64;
        REQUIRE( filter.bit_size() == bin_size_bits * optimal_bins );
    }
}

void validate_elements( std::string           output_filter_file,
                        const std::uint8_t    kmer_size,
                        const sequences_type& seqs,
                        const bins_type&      bins )
{
    // check if elements were properly inserted in the IBF
    // expects unique k-mers among all sequences without errors
    seqan3::interleaved_bloom_filter<> filter       = aux::load_ibf( output_filter_file );
    auto                               hash_adaptor = seqan3::views::kmer_hash( seqan3::ungapped{ kmer_size } );
    auto                               agent        = filter.counting_agent();

    std::vector< uint16_t >             expected_output( filter.bin_count(), 0 );
    seqan3::counting_vector< uint16_t > output( filter.bin_count(), 0 );

    int i = 0;
    for ( auto& seq : seqs )
    {
        // query IBF
        output += agent.bulk_count( seq | hash_adaptor );
        // Calculate expected (min) number of subsequences to be found (no errors)
        expected_output[bins[i]] = std::ranges::size( seq ) - kmer_size + 1;
        i += 1;
    }
    REQUIRE( output == expected_output );
}

void write_fasta( const std::string file, const sequences_type& seqs, const ids_type& ids )
{
    seqan3::sequence_file_output fout{ file };
    int                          i = 0;
    for ( auto& seq : seqs )
    {
        sequence_record_type rec{ seq, ids[i] };
        fout.push_back( rec );
        i += 1;
    }
}

void write_seqid_bin( std::string file, const sequences_type& seqs, const ids_type& ids, const bins_type& bins )
{
    // generate basic seqid_bin -> every sequence in one bin, no fragmentation
    std::ofstream seqid_bin_file{ file };
    uint16_t      i = 0;
    for ( auto& seq : seqs )
    {
        seqid_bin_file << ids[i] << "\t1\t" << std::ranges::size( seq ) << "\t" << bins[i] << '\n';
        i += 1;
    }
    seqid_bin_file.close();
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
        write_fasta( prefix + ".fasta", seqs, ids );
        cfg.reference_files = { prefix + ".fasta" };
        write_seqid_bin( prefix + "_seqid_bin.tsv", seqs, ids, bins );
        cfg.seqid_bin_file = prefix + "_seqid_bin.tsv";
    }
    return cfg;
}

} // namespace config_build


SCENARIO( "building indices", "[ganon-build]" )
{

    SECTION( "with default conf." )
    {
        auto cfg = config_build::defaultConfig( "default", seqs, ids, bins );
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, seqs, bins );
    }

    SECTION( "with multiple --reference-files" )
    {
        std::string prefix = "reference_files";
        auto        cfg    = config_build::defaultConfig( prefix );

        // write files separetly
        config_build::write_fasta( prefix + ".A.fasta", seqs, ids );
        config_build::write_fasta( prefix + ".B.fasta", extra_seqs, extra_ids );

        // merge entries
        auto merged_ids  = aux::vconcat( ids, extra_ids );
        auto merged_seqs = aux::vconcat( seqs, extra_seqs );
        auto merged_bins = aux::vconcat( bins, extra_bins );
        config_build::write_seqid_bin( prefix + "_seqid_bin.tsv", merged_seqs, merged_ids, merged_bins );

        cfg.reference_files = { prefix + ".A.fasta", prefix + ".B.fasta" };
        cfg.seqid_bin_file  = prefix + "_seqid_bin.tsv";

        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, merged_bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, merged_seqs, merged_bins );
    }

    SECTION( "with multiple --reference-files without --seqid-bin-file" )
    {
        // without --seqid-bin-file it should create one bin per sequence file
        std::string prefix = "reference_files_wo_seqid_bin";
        auto        cfg    = config_build::defaultConfig( prefix );

        // write one sequence per file
        config_build::write_fasta( prefix + ".S1.fasta", { seqs[0] }, { ids[0] } );
        config_build::write_fasta( prefix + ".S2.fasta", { seqs[1] }, { ids[1] } );
        config_build::write_fasta( prefix + ".S3.fasta", { seqs[2] }, { ids[2] } );
        cfg.reference_files = { prefix + ".S1.fasta", prefix + ".S2.fasta", prefix + ".S3.fasta" };

        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, seqs, bins );
    }

    SECTION( "with --kmer-size 11" )
    {
        auto cfg      = config_build::defaultConfig( "kmer_size_11", seqs, ids, bins );
        cfg.kmer_size = 11;
        // run ganon-build
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, seqs, bins );
    }

    SECTION( "with --kmer-size 27" )
    {
        auto cfg      = config_build::defaultConfig( "kmer_size_27", seqs, ids, bins );
        cfg.kmer_size = 27;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, seqs, bins );
    }

    SECTION( "with --kmer-size 31" )
    {
        // Should failed due to size limitation with dna5 (max 27, otherwise 32 with dna4)
        // https://docs.seqan.de/seqan/3-master-user/group__search__views.html#ga6e598d6a021868f704d39df73252974f
        auto cfg      = config_build::defaultConfig( "kmer_size_31", seqs, ids, bins );
        cfg.kmer_size = 31;
        REQUIRE_THROWS( GanonBuild::run( cfg ) );
    }

    SECTION( "with --filter-size-mb 2" )
    {
        auto cfg           = config_build::defaultConfig( "filter_size_mb_2", seqs, ids, bins );
        cfg.bin_size_bits  = 0;
        cfg.filter_size_mb = 2;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, seqs, bins );
    }

    SECTION( "with --hash-functions 2" )
    {
        auto cfg           = config_build::defaultConfig( "hash_functions_2", seqs, ids, bins );
        cfg.hash_functions = 2;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, seqs, bins );
    }

    SECTION( "with --directory-reference-files and --extension" )
    {
        // write file with specific extension "TEST.fasta"
        auto cfg            = config_build::defaultConfig( "directory_reference_files.TEST", seqs, ids, bins );
        cfg.reference_files = {};
        cfg.directory_reference_files = std::filesystem::canonical( "." );
        cfg.extension                 = ".TEST.fasta";
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter(
            cfg.output_filter_file, cfg.hash_functions, cfg.filter_size_mb, cfg.bin_size_bits, bins );
        config_build::validate_elements( cfg.output_filter_file, cfg.kmer_size, seqs, bins );
    }
}


SCENARIO( "updating indices", "[ganon-build]" )
{

    // build default filter
    auto cfg_build = config_build::defaultConfig( "update_base_build", seqs, ids, bins );
    REQUIRE( GanonBuild::run( cfg_build ) );
    config_build::validate_filter( cfg_build.output_filter_file,
                                   cfg_build.hash_functions,
                                   cfg_build.filter_size_mb,
                                   cfg_build.bin_size_bits,
                                   bins );
    config_build::validate_elements( cfg_build.output_filter_file, cfg_build.kmer_size, seqs, bins );

    SECTION( "with --update-filter-file creating 3 new bins" )
    {
        // update it
        auto cfg_update = config_build::defaultConfig( "3new_update", extra_seqs, extra_ids, extra_bins );
        // set filter to update
        cfg_update.update_filter_file = cfg_build.output_filter_file;
        REQUIRE( GanonBuild::run( cfg_update ) );

        auto merged_seqs = aux::vconcat( seqs, extra_seqs );
        auto merged_bins = aux::vconcat( bins, extra_bins );

        config_build::validate_filter( cfg_update.output_filter_file,
                                       cfg_update.hash_functions,
                                       cfg_update.filter_size_mb,
                                       cfg_update.bin_size_bits,
                                       merged_bins );
        config_build::validate_elements(
            cfg_update.output_filter_file, cfg_update.kmer_size, merged_seqs, merged_bins );

        SECTION( "with --update-complete" )
        {
            auto merged_ids = aux::vconcat( ids, extra_ids );
            // update complete, send all sequences
            auto cfg_update_complete =
                config_build::defaultConfig( "3new_update_complete", merged_seqs, merged_ids, merged_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;
            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete.output_filter_file,
                                           cfg_update_complete.hash_functions,
                                           cfg_update_complete.filter_size_mb,
                                           cfg_update_complete.bin_size_bits,
                                           merged_bins );
            config_build::validate_elements(
                cfg_update_complete.output_filter_file, cfg_update_complete.kmer_size, merged_seqs, merged_bins );

            // should be the same as not complete
            REQUIRE( aux::filesAreEqual( cfg_update_complete.output_filter_file, cfg_update.output_filter_file ) );
        }
    }

    SECTION( "with --update-filter-file creating 80 new bins" )
    {
        // set sequences to new bins
        const bins_type new_bins{ 12, 64, 82 };
        // update it
        auto cfg_update = config_build::defaultConfig( "80new_update", extra_seqs, extra_ids, new_bins );
        // set filter to update
        cfg_update.update_filter_file = cfg_build.output_filter_file;
        REQUIRE( GanonBuild::run( cfg_update ) );

        auto merged_seqs = aux::vconcat( seqs, extra_seqs );
        auto merged_bins = aux::vconcat( bins, new_bins );

        config_build::validate_filter( cfg_update.output_filter_file,
                                       cfg_update.hash_functions,
                                       cfg_update.filter_size_mb,
                                       cfg_update.bin_size_bits,
                                       merged_bins );
        config_build::validate_elements(
            cfg_update.output_filter_file, cfg_update.kmer_size, merged_seqs, merged_bins );

        // new file should be bigger (double bits + overhead)
        REQUIRE( aux::fileSizeBytes( cfg_update.output_filter_file )
                 > aux::fileSizeBytes( cfg_build.output_filter_file ) );

        SECTION( "with --update-complete" )
        {
            auto merged_ids = aux::vconcat( ids, extra_ids );
            // update complete, send all sequences
            auto cfg_update_complete =
                config_build::defaultConfig( "80new_update_complete", merged_seqs, merged_ids, merged_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;
            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete.output_filter_file,
                                           cfg_update_complete.hash_functions,
                                           cfg_update_complete.filter_size_mb,
                                           cfg_update_complete.bin_size_bits,
                                           merged_bins );
            config_build::validate_elements(
                cfg_update_complete.output_filter_file, cfg_update_complete.kmer_size, merged_seqs, merged_bins );

            // should be the same as not complete
            REQUIRE( aux::filesAreEqual( cfg_update_complete.output_filter_file, cfg_update.output_filter_file ) );
        }
    }

    SECTION( "with --update-filter-file without adding new bins (updating current bins)" )
    {
        // update it
        // use bins instead of extra_bins
        auto cfg_update = config_build::defaultConfig( "0new_update", extra_seqs, extra_ids, bins );
        // set filter to update
        cfg_update.update_filter_file = cfg_build.output_filter_file;
        REQUIRE( GanonBuild::run( cfg_update ) );
        config_build::validate_filter( cfg_update.output_filter_file,
                                       cfg_update.hash_functions,
                                       cfg_update.filter_size_mb,
                                       cfg_update.bin_size_bits,
                                       bins );
        // validate seqs
        config_build::validate_elements( cfg_update.output_filter_file, cfg_update.kmer_size, seqs, bins );
        // validate extra_seqs
        config_build::validate_elements( cfg_update.output_filter_file, cfg_update.kmer_size, extra_seqs, bins );

        SECTION( "with --update-complete" )
        {
            auto merged_seqs   = aux::vconcat( seqs, extra_seqs );
            auto merged_ids    = aux::vconcat( ids, extra_ids );
            auto repeated_bins = aux::vconcat( bins, bins );
            // update complete, send all sequences
            auto cfg_update_complete =
                config_build::defaultConfig( "0new_update_complete", merged_seqs, merged_ids, repeated_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;
            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete.output_filter_file,
                                           cfg_update_complete.hash_functions,
                                           cfg_update_complete.filter_size_mb,
                                           cfg_update_complete.bin_size_bits,
                                           bins );

            // validate seqs
            config_build::validate_elements(
                cfg_update_complete.output_filter_file, cfg_update_complete.kmer_size, seqs, bins );
            // validate extra_seqs
            config_build::validate_elements(
                cfg_update_complete.output_filter_file, cfg_update_complete.kmer_size, extra_seqs, bins );

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

            auto cfg_update_complete =
                config_build::defaultConfig( "0new_update_complete_remove", merged_seqs, merged_ids, merged_bins );
            cfg_update_complete.update_filter_file = cfg_build.output_filter_file;
            cfg_update_complete.update_complete    = true;

            REQUIRE( GanonBuild::run( cfg_update_complete ) );
            config_build::validate_filter( cfg_update_complete.output_filter_file,
                                           cfg_update_complete.hash_functions,
                                           cfg_update_complete.filter_size_mb,
                                           cfg_update_complete.bin_size_bits,
                                           merged_bins );

            // validate seqs (without first)
            config_build::validate_elements(
                cfg_update_complete.output_filter_file, cfg_update_complete.kmer_size, rem_seqs, rem_bins );
            // validate extra_seqs (all should be there)
            config_build::validate_elements(
                cfg_update_complete.output_filter_file, cfg_update_complete.kmer_size, extra_seqs, bins );

            // should be different as before
            REQUIRE_FALSE(
                aux::filesAreEqual( cfg_update_complete.output_filter_file, cfg_update.output_filter_file ) );
        }
    }
}
