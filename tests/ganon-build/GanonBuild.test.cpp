#include "aux/Aux.hpp"

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

namespace config_build
{

GanonBuild::Config defaultConfig()
{
    GanonBuild::Config cfg;
    cfg.seqid_bin_file     = "filters/bacteria_acc_bin.txt";
    cfg.output_filter_file = "test_output.ibf";
    cfg.filter_size_bits   = 8388352;
    cfg.reference_files    = { "sequences/bacteria_NC_010333.1.fasta.gz",
                            "sequences/bacteria_NC_017163.1.fasta.gz",
                            "sequences/bacteria_NC_017164.1.fasta.gz",
                            "sequences/bacteria_NC_017543.1.fasta.gz" };
    cfg.verbose            = false;
    cfg.quiet              = true;
    return cfg;
}

const std::string bacteria_filter = "filters/bacteria_build.ibf";

} // namespace config_build

SCENARIO( "Build", "[ganon-build]" )
{
    auto cfg = config_build::defaultConfig();

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_filter_file, config_build::bacteria_filter ) );
}

SCENARIO( "Build from folder", "[ganon-build]" )
{
    auto cfg                      = config_build::defaultConfig();
    cfg.directory_reference_files = "sequences/";
    cfg.extension                 = ".gz";

    REQUIRE( GanonBuild::run( cfg ) );
}

SCENARIO( "Build forced failure with different k-mer size", "[ganon-build]" )
{
    auto cfg      = config_build::defaultConfig();
    cfg.kmer_size = 18;

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::bacteria_filter ) );
}

SCENARIO( "Build forced failure with different number of hash functions", "[ganon-build]" )
{
    auto cfg           = config_build::defaultConfig();
    cfg.hash_functions = 2;

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::bacteria_filter ) );
}

SCENARIO( "Build forced failure with different bloom size", "[ganon-build]" )
{
    auto cfg             = config_build::defaultConfig();
    cfg.filter_size_bits = 17000000;

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::bacteria_filter ) );
}

SCENARIO( "Build forced failure with different reference files", "[ganon-build]" )
{
    auto cfg = config_build::defaultConfig();
    cfg.reference_files.erase( cfg.reference_files.begin() );

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::bacteria_filter ) );
}
SCENARIO( "Build forced failure with incomplete seqid-bin file", "[ganon-build]" )
{
    auto cfg           = config_build::defaultConfig();
    cfg.seqid_bin_file = "filters/bacteria_acc_bin_incomplete.txt";

    REQUIRE( GanonBuild::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_filter_file, config_build::bacteria_filter ) );
}

SCENARIO( "Update", "[ganon-build]" )
{
    // update bacteria filter with virus
    auto cfg               = config_build::defaultConfig();
    cfg.seqid_bin_file     = "filters/bacteria_upd_virus_acc_bin.txt";
    cfg.reference_files    = { "sequences/virus_NC_003676.1.fasta.gz",
                            "sequences/virus_NC_011646.1.fasta.gz",
                            "sequences/virus_NC_032412.1.fasta.gz",
                            "sequences/virus_NC_035470.1.fasta.gz" };
    cfg.update_filter_file = config_build::bacteria_filter;

    REQUIRE( GanonBuild::run( cfg ) );

    // should have the same size since number of new bins does not pass the multiple of 64 threshold
    // bacteria 67 bins + virus 9 bins
    REQUIRE( aux::fileSize( cfg.output_filter_file ) == aux::fileSize( cfg.update_filter_file ) );

    // check file
    REQUIRE( aux::filesAreEqual( cfg.output_filter_file, "filters/bacteria_upd_virus_build.ibf" ) );
}

SCENARIO( "Update complete", "[ganon-build]" )
{
    // update complete bacteria filter with virus
    // removing NC_017164.1 bin 63
    auto cfg               = config_build::defaultConfig();
    cfg.seqid_bin_file     = "filters/bacteria_upd_complete_virus_acc_bin.txt";
    cfg.reference_files    = { "sequences/bacteria_NC_010333.1.fasta.gz", "sequences/bacteria_NC_017163.1.fasta.gz",
                            "sequences/bacteria_NC_017543.1.fasta.gz", "sequences/virus_NC_003676.1.fasta.gz",
                            "sequences/virus_NC_011646.1.fasta.gz",    "sequences/virus_NC_032412.1.fasta.gz",
                            "sequences/virus_NC_035470.1.fasta.gz" };
    cfg.update_complete    = true;
    cfg.update_filter_file = config_build::bacteria_filter;

    REQUIRE( GanonBuild::run( cfg ) );

    // should have the same size since number of new bins does not pass the multiple of 64 threshold
    // bacteria 67 bins + virus 8 bins
    REQUIRE( aux::fileSize( cfg.output_filter_file ) == aux::fileSize( cfg.update_filter_file ) );

    // check file
    REQUIRE( aux::filesAreEqual( cfg.output_filter_file, "filters/bacteria_upd_complete_virus_build.ibf" ) );
}

SCENARIO( "Update increasing filter size", "[ganon-build]" )
{
    // update virus filter with bacteria and archaea
    auto cfg            = config_build::defaultConfig();
    cfg.seqid_bin_file  = "filters/virus_upd_archaea_bacteria_acc_bin.txt";
    cfg.reference_files = {
        "sequences/archaea_NC_015430.1.fasta.gz",  "sequences/archaea_NC_023011.1.fasta.gz",
        "sequences/archaea_NC_023012.1.fasta.gz",  "sequences/bacteria_NC_010333.1.fasta.gz",
        "sequences/bacteria_NC_017163.1.fasta.gz", "sequences/bacteria_NC_017164.1.fasta.gz",
        "sequences/bacteria_NC_017543.1.fasta.gz",
    };

    cfg.update_filter_file = "filters/virus.ibf";

    REQUIRE( GanonBuild::run( cfg ) );

    // should have increased the filter size to keep the original FP
    REQUIRE( aux::fileSize( cfg.output_filter_file ) > aux::fileSize( cfg.update_filter_file ) );

    // check file
    REQUIRE( aux::filesAreEqual( cfg.output_filter_file, "filters/virus_upd_archaea_bacteria.ibf" ) );
}