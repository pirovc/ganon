#include "aux/Aux.hpp"

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

namespace detail
{

Config testConfig()
{
    Config cfg;
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

    return cfg;
}

const std::string outputFile = "build_output.filter";

} // namespace detail

SCENARIO( "Bacteria", "[ganon-build]" )
{
    const auto cfg = detail::testConfig();

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_filter_file, detail::outputFile ) );
}

SCENARIO( "Bacteria forced failure", "[ganon-build]" )
{
    auto cfg      = detail::testConfig();
    cfg.kmer_size = 18;

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, detail::outputFile ) );
}
