#include "aux/Aux.hpp"

#include <ganon-classify/Config.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

namespace config_classify
{

GanonClassify::Config defaultConfig()
{
    GanonClassify::Config cfg;
    cfg.output_single       = true;
    cfg.output_all          = true;
    cfg.output_unclassified = false;
    cfg.kmer_size           = { 19 };
    cfg.threads             = 4;
    cfg.verbose             = false;
    cfg.quiet               = true;
    cfg.hierarchy_labels    = { "1" };
    cfg.offset              = 1;

    return cfg;
}

std::vector< std::string > output_ext{ "all", "lca", "rep" };
std::string                results_path = "results/";

} // namespace config_classify

// Static results

SCENARIO( "Classify", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error     = { 3 };
    cfg.output_prefix = "b-b_e3";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify LCA", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error     = { 0 }; // one match per read - all=lca
    cfg.output_prefix = "b-b_e0";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + ".all", cfg.output_prefix + ".lca" ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }

    cfg.max_error     = { 5 }; // more than one match per read
    cfg.output_prefix = "b-b_e5";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqualSorted( cfg.output_prefix + ".all", cfg.output_prefix + ".lca" ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify paired-reads concat", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.paired_reads  = { "reads/bacteria_id.1.fq", "reads/bacteria_id.2.fq" };
    cfg.max_error     = { 0 };
    cfg.output_prefix = "b-b_e0-paired";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify paired-reads concat with unique errors", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf" };
    cfg.map              = { "filters/bacteria.map" };
    cfg.tax              = { "filters/bacteria.tax" };
    cfg.paired_reads     = { "reads/bacteria_id.1.fq", "reads/bacteria_id.2.fq" };
    cfg.max_error        = { 1 };
    cfg.max_error_unique = { 0 };
    cfg.output_prefix    = "b-b_e1u0-paired";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify reads with no filtering", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error     = { 5 };
    cfg.strata_filter = { -1 };
    cfg.output_prefix = "b-b_e5l-1";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify paired-reads and single-reads with multiple indices", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/archaea.ibf", "filters/bacteria.ibf" };
    cfg.map              = { "filters/archaea.map", "filters/bacteria.map" };
    cfg.tax              = { "filters/archaea.tax", "filters/bacteria.tax" };
    cfg.paired_reads     = { "reads/bacteria_id.1.fq", "reads/bacteria_id.2.fq" };
    cfg.single_reads     = { "reads/archaea.simulated.1.fq" };
    cfg.max_error        = { 0 };
    cfg.hierarchy_labels = { "1", "2" };
    cfg.output_prefix    = "ab-ab_e0-paired";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

/*SCENARIO( "Classify with offset", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.offset        = 2;
    cfg.max_error     = { 3 };
    cfg.output_prefix = "b-b_e3f2";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}*/

SCENARIO( "Classify with no errors allowed", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error     = { 0 };
    cfg.output_prefix = "b-b_e0-noerrors";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }

    // using min_kmers = 1 should achieve the same result
    auto cfg2          = config_classify::defaultConfig();
    cfg2.ibf           = { "filters/bacteria.ibf" };
    cfg2.map           = { "filters/bacteria.map" };
    cfg2.tax           = { "filters/bacteria.tax" };
    cfg2.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg2.min_kmers     = { 1 };
    cfg2.output_prefix = "b-b_m1-noerrors";

    INFO( "output_prefix2: " + cfg2.output_prefix );
    REQUIRE( GanonClassify::run( cfg2 ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg2.output_prefix + "." + ext,
                                           config_classify::results_path + cfg2.output_prefix + "." + ext ) );
    }

    // check if both files are the same
    for ( auto const& ext : config_classify::output_ext )
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext, cfg2.output_prefix + "." + ext ) );
}

SCENARIO( "Classify with min kmers", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.min_kmers     = { 0.29 }; // should work the same as -e 3, threshold = 25 19-mers
    cfg.output_prefix = "b-b_m0.29";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify with different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf" };
    cfg.map              = { "filters/bacteria.map" };
    cfg.tax              = { "filters/bacteria.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error_unique = { 1 };
    cfg.max_error        = { 3 };
    cfg.output_prefix    = "b-b_e3u1";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

/*SCENARIO( "Classify with offset and different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf" };
    cfg.map              = { "filters/bacteria.map" };
    cfg.tax              = { "filters/bacteria.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error        = { 2 };
    cfg.max_error_unique = { 0 };
    cfg.offset           = 6;
    cfg.output_prefix    = "b-b_e2u0f6";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}*/

/*SCENARIO( "Classify with offset, min kmers and different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf" };
    cfg.map              = { "filters/bacteria.map" };
    cfg.tax              = { "filters/bacteria.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq" };
    cfg.min_kmers        = { 0.53 }; // should be equal as 2 errors (9 k-mers)
    cfg.max_error_unique = { 0 };
    cfg.offset           = 6;
    cfg.output_prefix    = "b-b_m0.53u0f6";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}*/

/*SCENARIO( "Classify with offset higher than k", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.offset        = 25; // should be limited by k-mer size (19)
    cfg.output_prefix = "b-b_f25";
    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );

    auto cfg2          = config_classify::defaultConfig();
    cfg2.ibf           = { "filters/bacteria.ibf" };
    cfg2.map           = { "filters/bacteria.map" };
    cfg2.tax           = { "filters/bacteria.tax" };
    cfg2.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg2.offset        = 19;
    cfg2.output_prefix = "b-b_f19";
    INFO( "output_prefix2: " + cfg2.output_prefix );
    REQUIRE( GanonClassify::run( cfg2 ) );

    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext, cfg2.output_prefix + "." + ext ) );
    }
}*/

SCENARIO( "Classify multi-filter without errors allowed", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map           = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax           = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error     = { 0 };
    cfg.output_prefix = "ba-ba_e0";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-filter with errors allowed", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map           = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax           = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error     = { 4 };
    cfg.output_prefix = "ba-ba_e4";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-filter with multiple errors", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map           = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax           = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error     = { 0, 4 };
    cfg.output_prefix = "ba-ba_e04";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-hierarchy without errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map              = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax              = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error        = { 0 };
    cfg.hierarchy_labels = { "1", "2" };
    cfg.output_prefix    = "ba-ba_e0c12";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-hierarchy split files", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map              = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax              = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error        = { 0 };
    cfg.hierarchy_labels = { "1", "2" };
    cfg.output_prefix    = "ba-ba_e0c12-split";
    std::string output1  = cfg.output_prefix + ".1.all";
    std::string output2  = cfg.output_prefix + ".2.all";
    cfg.output_single    = false;

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + ".rep",
                                       config_classify::results_path + cfg.output_prefix + ".rep" ) );
    REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + ".1.lca",
                                       config_classify::results_path + cfg.output_prefix + ".1.lca" ) );
    REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + ".2.lca",
                                       config_classify::results_path + cfg.output_prefix + ".2.lca" ) );
    REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + ".1.all",
                                       config_classify::results_path + cfg.output_prefix + ".1.all" ) );
    REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + ".2.all",
                                       config_classify::results_path + cfg.output_prefix + ".2.all" ) );
}

SCENARIO( "Classify multi-hierarchy with errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map              = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax              = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error        = { 4 };
    cfg.hierarchy_labels = { "1", "2" };
    cfg.output_prefix    = "ba-ba_e4c12";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-hierarchy with multiple errors", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/archaea.ibf", "filters/bacteria.ibf" };
    cfg.map              = { "filters/archaea.map", "filters/bacteria.map" };
    cfg.tax              = { "filters/archaea.tax", "filters/bacteria.tax" };
    cfg.single_reads     = { "reads/archaea.simulated.1.fq" };
    cfg.hierarchy_labels = { "1", "2" };
    cfg.max_error        = { 3, 4 };
    cfg.output_prefix    = "ab-a_e34c12";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-hierarchy with multiple errors and multiple unique errors", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/archaea.ibf", "filters/bacteria.ibf" };
    cfg.map              = { "filters/archaea.map", "filters/bacteria.map" };
    cfg.tax              = { "filters/archaea.tax", "filters/bacteria.tax" };
    cfg.single_reads     = { "reads/archaea.simulated.1.fq" };
    cfg.hierarchy_labels = { "1", "2" };
    cfg.max_error        = { 3, 4 };
    cfg.max_error_unique = { 0, 1 };
    cfg.output_prefix    = "ab-a_e34c12u01";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}

// Functionality

/*SCENARIO( "Classify problematic fastq", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/problematic.fq" };
    cfg.max_error     = { 3 };
    cfg.output_prefix = "b-problematic_e3";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileLines( cfg.output_prefix + ".all" ) == 4 );

    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg.output_prefix + "." + ext,
                                           config_classify::results_path + cfg.output_prefix + "." + ext ) );
    }
}*/

SCENARIO( "Classify without matches", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/virus.simulated.1.fq" };
    cfg.max_error     = { 0 };
    cfg.output_prefix = "b-v_e0";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );
}

SCENARIO( "Classify multi-filter without matches", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map           = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax           = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads  = { "reads/virus.simulated.1.fq" };
    cfg.max_error     = { 0 };
    cfg.output_prefix = "ba-v_e0";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );
}

SCENARIO( "Classify multi-hierarchy without matches", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.ibf", "filters/archaea.ibf" };
    cfg.map              = { "filters/bacteria.map", "filters/archaea.map" };
    cfg.tax              = { "filters/bacteria.tax", "filters/archaea.tax" };
    cfg.single_reads     = { "reads/virus.simulated.1.fq" };
    cfg.max_error        = { 0 };
    cfg.hierarchy_labels = { "1", "2" };
    cfg.output_prefix    = "ba-v_e0c12";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );
}

SCENARIO( "Classify forced failure with different max. errors allowed", "[ganon-classify]" )
{
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error     = { 0 };
    cfg.output_prefix = "b-b_e0-failure";

    const std::string undesired_output_prefix = "b-b_e3";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE_FALSE( aux::filesAreEqualSorted(
            cfg.output_prefix + "." + ext, config_classify::results_path + undesired_output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-filter with partial matching reads", "[ganon-classify]" )
{
    // should match results as if working with only one filter because max_error = 0
    auto cfg1          = config_classify::defaultConfig();
    cfg1.ibf           = { "filters/archaea.ibf" };
    cfg1.map           = { "filters/archaea.map" };
    cfg1.tax           = { "filters/archaea.tax" };
    cfg1.single_reads  = { "reads/archaea.simulated.1.fq" };
    cfg1.max_error     = { 0 };
    cfg1.output_prefix = "a-a_e0";

    auto cfg2          = config_classify::defaultConfig();
    cfg2.ibf           = { "filters/bacteria.ibf", "filters/archaea.ibf", "filters/virus.ibf" };
    cfg2.map           = { "filters/bacteria.map", "filters/archaea.map", "filters/virus.map" };
    cfg2.tax           = { "filters/bacteria.tax", "filters/archaea.tax", "filters/virus.tax" };
    cfg2.single_reads  = { "reads/archaea.simulated.1.fq" };
    cfg2.max_error     = { 0 };
    cfg2.output_prefix = "bav-a_e0";

    INFO( "output_prefix: " + cfg1.output_prefix );
    INFO( "output_prefix2: " + cfg2.output_prefix );
    REQUIRE( GanonClassify::run( cfg1 ) );
    REQUIRE( GanonClassify::run( cfg2 ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg1.output_prefix + "." + ext, cfg2.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify multi-hierarchy with partial matching reads", "[ganon-classify]" )
{
    // should match results as if working with only one filter because max_error = 0
    // reads should 'survive' till the last filter
    auto cfg1             = config_classify::defaultConfig();
    cfg1.ibf              = { "filters/virus.ibf" };
    cfg1.map              = { "filters/virus.map" };
    cfg1.tax              = { "filters/virus.tax" };
    cfg1.single_reads     = { "reads/virus.simulated.1.fq" };
    cfg1.max_error        = { 0 };
    cfg1.hierarchy_labels = { "3" }; // to match .rep
    cfg1.output_prefix    = "v-v_e0";

    // Test if the reads are surviving the hierachies
    auto cfg2             = config_classify::defaultConfig();
    cfg2.ibf              = { "filters/bacteria.ibf", "filters/archaea.ibf", "filters/virus.ibf" };
    cfg2.map              = { "filters/bacteria.map", "filters/archaea.map", "filters/virus.map" };
    cfg2.tax              = { "filters/bacteria.tax", "filters/archaea.tax", "filters/virus.tax" };
    cfg2.single_reads     = { "reads/virus.simulated.1.fq" };
    cfg2.max_error        = { 0 };
    cfg2.hierarchy_labels = { "1", "2", "3" };
    cfg2.output_prefix    = "bav-v_e0";

    INFO( "output_prefix: " + cfg1.output_prefix );
    REQUIRE( GanonClassify::run( cfg1 ) );
    INFO( "output_prefix2: " + cfg2.output_prefix );
    REQUIRE( GanonClassify::run( cfg2 ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg1.output_prefix + "." + ext, cfg2.output_prefix + "." + ext ) );
    }
}

SCENARIO( "Classify after update", "[ganon-classify]" )
{
    // No match
    auto cfg          = config_classify::defaultConfig();
    cfg.ibf           = { "filters/bacteria.ibf" };
    cfg.map           = { "filters/bacteria.map" };
    cfg.tax           = { "filters/bacteria.tax" };
    cfg.single_reads  = { "reads/virus.simulated.1.fq" };
    cfg.max_error     = { 3 };
    cfg.output_prefix = "b-v_e3";

    INFO( "output_prefix: " + cfg.output_prefix );
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );

    GanonBuild::Config cfg_build;
    cfg_build.update_filter_file = "filters/bacteria.ibf";
    cfg_build.seqid_bin_file     = "filters/bacteria_upd_virus_acc_bin.txt";
    cfg_build.output_filter_file = "bacteria_virus.ibf";
    // cfg_build.filter_size_bits   = 8388608;
    cfg_build.bin_size_bits   = 65534;
    cfg_build.reference_files = { "sequences/virus_NC_003676.1.fasta.gz",
                                  "sequences/virus_NC_011646.1.fasta.gz",
                                  "sequences/virus_NC_032412.1.fasta.gz",
                                  "sequences/virus_NC_035470.1.fasta.gz" };
    cfg_build.verbose         = cfg.verbose;
    cfg_build.quiet           = cfg.quiet;

    REQUIRE( GanonBuild::run( cfg_build ) );

    auto cfg2                               = config_classify::defaultConfig();
    cfg2.max_error                          = { 3 };
    cfg2.ibf                                = { cfg_build.output_filter_file };
    cfg2.map                                = { "filters/bacteria_virus.map" };
    cfg2.tax                                = { "filters/bacteria_virus.tax" };
    cfg2.single_reads                       = { "reads/bacteria.simulated.1.fq" };
    cfg2.output_prefix                      = "bv-b_e3";
    const std::string desired_output_prefix = "b-b_e3";

    // Results of bacteria should be the same as without update (check if FP changed or filter was affected)
    INFO( "output_prefix2: " + cfg2.output_prefix );
    REQUIRE( GanonClassify::run( cfg2 ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg2.output_prefix + "." + ext,
                                           config_classify::results_path + desired_output_prefix + "." + ext ) );
    }

    // Virus reads should now map
    auto cfg3          = config_classify::defaultConfig();
    cfg3.max_error     = { 3 };
    cfg3.ibf           = { cfg_build.output_filter_file };
    cfg3.map           = { "filters/bacteria_virus.map" };
    cfg3.tax           = { "filters/bacteria_virus.tax" };
    cfg3.single_reads  = { "reads/virus.simulated.1.fq" };
    cfg3.output_prefix = "bv-v_e3";

    INFO( "output_prefix2: " + cfg3.output_prefix );
    REQUIRE( GanonClassify::run( cfg3 ) );
    for ( auto const& ext : config_classify::output_ext )
    {
        INFO( "extension: " + ext );
        REQUIRE( aux::filesAreEqualSorted( cfg3.output_prefix + "." + ext,
                                           config_classify::results_path + cfg3.output_prefix + "." + ext ) );
    }
}
