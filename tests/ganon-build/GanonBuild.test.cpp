#include "aux/Aux.hpp"

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

namespace config_build
{

GanonBuild::Config defaultConfig()
{
    GanonBuild::Config cfg;
    cfg.seqid_bin_file     = "bacteria_seqid_bin.txt";
    cfg.output_filter_file = "test_output.filter";
    cfg.filter_size        = 15797760;
    cfg.kmer_size          = 19;
    cfg.hash_functions     = 3;
    cfg.reference_files    = { "bacteria_NC_010333.1.fasta.gz",
                            "bacteria_NC_017163.1.fasta.gz",
                            "bacteria_NC_017164.1.fasta.gz",
                            "bacteria_NC_017543.1.fasta.gz" };
    cfg.threads            = 1;
    cfg.build_threads      = 1;
    cfg.verbose            = true;

    return cfg;
}

const std::string outputFile = "build_output.filter";

} // namespace config_build

SCENARIO( "Build", "[ganon-build]" )
{
    const auto cfg = config_build::defaultConfig();

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_filter_file, config_build::outputFile ) );
}

SCENARIO( "Build forced failure with different k-mer size", "[ganon-build]" )
{
    auto cfg      = config_build::defaultConfig();
    cfg.kmer_size = 18;

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::outputFile ) );
}

SCENARIO( "Build forced failure with different number of hash functions", "[ganon-build]" )
{
    auto cfg           = config_build::defaultConfig();
    cfg.hash_functions = 2;

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::outputFile ) );
}

SCENARIO( "Build forced failure with different bloom size", "[ganon-build]" )
{
    auto cfg        = config_build::defaultConfig();
    cfg.filter_size = 17000000;

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::outputFile ) );
}

SCENARIO( "Build forced failure with different reference files", "[ganon-build]" )
{
    auto cfg = config_build::defaultConfig();
    cfg.reference_files.erase( cfg.reference_files.begin() );

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::outputFile ) );
}
SCENARIO( "Build forced failure with incomplete seqid-bin file", "[ganon-build]" )
{
    auto cfg           = config_build::defaultConfig();
    cfg.seqid_bin_file = "bacteria_seqid_bin_incomplete.txt";

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::outputFile ) );
}