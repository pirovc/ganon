#include "aux/Aux.hpp"

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

SCENARIO( "Bacteria", "[ganon-build]" )
{
    Config config;
    config.seqid_bin_file     = "bacteria_seqid_bin.txt";
    config.output_filter_file = "test_output.filter";
    config.filter_size        = 15797760;
    config.kmer_size          = 19;
    config.hash_functions     = 3;
    config.reference_files    = { "bacteria_NC_010333.1.fasta.gz",
                               "bacteria_NC_017163.1.fasta.gz",
                               "bacteria_NC_017164.1.fasta.gz",
                               "bacteria_NC_017543.1.fasta.gz" };
    config.threads            = 1;
    config.build_threads      = 1;

    REQUIRE( GanonBuild::run( config ) );
    REQUIRE( aux::filesAreEqual( config.output_filter_file, "build_output.filter" ) );
}
