#include "aux/Aux.hpp"


#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>

#include <seqan3/std/ranges>

#include <filesystem>
#include <iostream>

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

using namespace seqan3::literals;

using sequences_type = std::vector< seqan3::dna5_vector >;

namespace config_build
{


using sequence_record_type = seqan3::sequence_record< seqan3::type_list< std::vector< seqan3::dna5 >, std::string >,
                                                      seqan3::fields< seqan3::field::seq, seqan3::field::id > >;


void validate_filter( const GanonBuild::Config& cfg, sequences_type& seqs )
{
    // validate generated filter file in properties and query

    // check if file exists
    REQUIRE( std::filesystem::exists( cfg.output_filter_file ) );
    // load filter
    seqan3::interleaved_bloom_filter<> filter = aux::load_ibf( cfg.output_filter_file );
    // check bin count
    REQUIRE( filter.bin_count() == seqs.size() );
    // check hash functions
    REQUIRE( filter.hash_function_count() == cfg.hash_functions );

    // check size
    if ( cfg.filter_size_mb > 0 )
    {
        REQUIRE( filter.bit_size() == cfg.filter_size_mb * 8388608u );
    }
    else
    {
        uint64_t optimal_bins = ( std::floor( filter.bin_count() / 64 ) + 1 ) * 64;
        REQUIRE( filter.bit_size() == cfg.bin_size_bits * optimal_bins );
    }

    // check elements inserted
    auto hash_adaptor = seqan3::views::kmer_hash( seqan3::ungapped{ cfg.kmer_size } );
    auto agent        = filter.counting_agent();

    std::vector< uint16_t >             expected_output( seqs.size(), 0 );
    seqan3::counting_vector< uint16_t > output( seqs.size(), 0 );

    int i = 0;
    for ( auto& seq : seqs )
    {
        // query IBF
        output += agent.bulk_count( seq | hash_adaptor );
        // Calculate expected (min) number of subsequences to be found (no errors)
        expected_output[i] = std::ranges::size( seq ) - cfg.kmer_size + 1;
        REQUIRE( output == expected_output );
        i += 1;
    }
}

void write_fasta( sequences_type& seqs, std::string file )
{
    seqan3::sequence_file_output fout{ file };
    int                          i = 0;
    for ( auto& seq : seqs )
    {
        sequence_record_type rec{ seq, "S" + std::to_string( i ) };
        fout.push_back( rec );
        i += 1;
    }
}

void write_seqid_bin( sequences_type& seqs, std::string file )
{
    // generate basic seqid_bin file with one full sequence for each bin
    std::ofstream seqid_bin_file{ file };
    uint16_t      i = 0;
    for ( auto& seq : seqs )
    {
        seqid_bin_file << "S" << i << "\t1\t" << std::ranges::size( seq ) << "\t" << i << '\n';
        i += 1;
    }
    seqid_bin_file.close();
}

GanonBuild::Config defaultConfig( std::string prefix, sequences_type& seqs )
{

    // Make config
    GanonBuild::Config cfg;
    cfg.bin_size_bits      = 5000;
    cfg.verbose            = false;
    cfg.quiet              = true;
    cfg.kmer_size          = 19;
    cfg.hash_functions     = 3;
    cfg.output_filter_file = prefix + ".ibf";

    // if input sequences are sent, create files
    if ( seqs.size() )
    {
        write_fasta( seqs, prefix + ".fasta" );
        cfg.reference_files = { prefix + ".fasta" };
        write_seqid_bin( seqs, prefix + "_seqid_bin.tsv" );
        cfg.seqid_bin_file = prefix + "_seqid_bin.tsv";
    }
    return cfg;
}

} // namespace config_build

SCENARIO( "building indices", "[ganon-build]" )
{

    // Sequences to build the ibf
    sequences_type seqs{ "TTCAATTCGGCGTACTCAGCATCGCAGCTAGCTGTACGGCTAGTCGTCAT"_dna5,
                         "TTGGGGCTAAACAGCACTATACAGGCGGCTAGCATGTATTAGGGGAGCTC"_dna5,
                         "ACCTTCGATTTCTTTAGATCGGGGATGATGATGCATGATGCTTAGGGATT"_dna5 };

    SECTION( "with default conf." )
    {
        auto cfg = config_build::defaultConfig( "default", seqs );
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, seqs );
    }


    SECTION( "with --kmer-size 11" )
    {
        auto cfg      = config_build::defaultConfig( "kmer_size_11", seqs );
        cfg.kmer_size = 11;
        // run ganon-build
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, seqs );
    }

    SECTION( "with --kmer-size 27" )
    {
        auto cfg      = config_build::defaultConfig( "kmer_size_27", seqs );
        cfg.kmer_size = 27;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, seqs );
    }

    SECTION( "with --kmer-size 31" )
    {
        // Should failed due to size limitation with dna5 (max 27, otherwise 32 with dna4)
        // https://docs.seqan.de/seqan/3-master-user/group__search__views.html#ga6e598d6a021868f704d39df73252974f
        auto cfg      = config_build::defaultConfig( "kmer_size_31", seqs );
        cfg.kmer_size = 31;
        REQUIRE_THROWS( GanonBuild::run( cfg ) );
    }

    SECTION( "with --filter-size-mb 2" )
    {
        auto cfg           = config_build::defaultConfig( "filter_size_mb_2", seqs );
        cfg.bin_size_bits  = 0;
        cfg.filter_size_mb = 2;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, seqs );
    }

    SECTION( "with --hash-functions 2" )
    {
        auto cfg           = config_build::defaultConfig( "hash_functions_2", seqs );
        cfg.hash_functions = 2;
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, seqs );
    }

    SECTION( "with --directory-reference-files and --extension" )
    {
        // write file with specific extension "TEST.fasta"
        auto cfg                      = config_build::defaultConfig( "directory_reference_files.TEST", seqs );
        cfg.reference_files           = {};
        cfg.directory_reference_files = std::filesystem::canonical( "." );
        cfg.extension                 = ".TEST.fasta";
        REQUIRE( GanonBuild::run( cfg ) );
        config_build::validate_filter( cfg, seqs );
    }
}
