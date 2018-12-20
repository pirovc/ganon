#include "aux/Aux.hpp"

#include <ganon-classify/Config.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <catch2/catch.hpp>

namespace config_classify
{

GanonClassify::Config defaultConfig()
{
    GanonClassify::Config cfg;
    cfg.output_file      = "classify_test_output.txt";
    cfg.max_error        = 3;
    cfg.max_error_unique = -1;
    cfg.verbose          = true;
    cfg.testing          = true;
    cfg.threads          = 3;
    return cfg;
}

} // namespace config_classify


// Static results

SCENARIO( "Classify", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter" };
    cfg.group_bin_files              = { "bacteria.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq" };
    const std::string desired_output = "classify_output-b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}

SCENARIO( "Classify with no errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter" };
    cfg.group_bin_files              = { "bacteria.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq" };
    cfg.max_error                    = 0;
    const std::string desired_output = "classify_output-b-b_e0.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}

SCENARIO( "Classify with different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg               = config_classify::defaultConfig();
    cfg.bloom_filter_files = { "bacteria.filter" };
    cfg.group_bin_files    = { "bacteria.map" };
    cfg.reads              = { "bacteria.simulated.1.fq" };
    cfg.max_error_unique   = 1;

    const std::string desired_output = "classify_output-b-b_e3u1.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}

SCENARIO( "Classify multi-filter without errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter", "archaea.filter" };
    cfg.group_bin_files              = { "bacteria.map", "archaea.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq", "archaea.simulated.1.fq" };
    cfg.max_error                    = 0;
    const std::string desired_output = "classify_output-ba-ba_e0.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}

SCENARIO( "Classify multi-filter with errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter", "archaea.filter" };
    cfg.group_bin_files              = { "bacteria.map", "archaea.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq", "archaea.simulated.1.fq" };
    cfg.max_error                    = 4;
    const std::string desired_output = "classify_output-ba-ba_e4.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}

SCENARIO( "Classify multi-hierarchy without errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter", "archaea.filter" };
    cfg.group_bin_files              = { "bacteria.map", "archaea.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq", "archaea.simulated.1.fq" };
    cfg.max_error                    = 0;
    cfg.filter_hierarchy             = "1,2";
    const std::string desired_output = "classify_output-ba-ba_e0c12.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}

SCENARIO( "Classify multi-hierarchy with errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter", "archaea.filter" };
    cfg.group_bin_files              = { "bacteria.map", "archaea.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq", "archaea.simulated.1.fq" };
    cfg.max_error                    = 4;
    cfg.filter_hierarchy             = "1,2";
    const std::string desired_output = "classify_output-ba-ba_e4c12.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}


// Functionality

SCENARIO( "Classify without matches", "[ganon-classify]" )
{
    auto cfg               = config_classify::defaultConfig();
    cfg.bloom_filter_files = { "bacteria.filter" };
    cfg.group_bin_files    = { "bacteria.map" };
    cfg.reads              = { "virus.simulated.1.fq" };
    cfg.max_error          = 0;
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_file ) );
}

SCENARIO( "Classify multi-filter without matches", "[ganon-classify]" )
{
    auto cfg               = config_classify::defaultConfig();
    cfg.bloom_filter_files = { "bacteria.filter", "archaea.filter" };
    cfg.group_bin_files    = { "bacteria.map", "archaea.map" };
    cfg.reads              = { "virus.simulated.1.fq" };
    cfg.max_error          = 0;
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_file ) );
}

SCENARIO( "Classify multi-hierarchy without matches", "[ganon-classify]" )
{
    auto cfg               = config_classify::defaultConfig();
    cfg.bloom_filter_files = { "bacteria.filter", "archaea.filter" };
    cfg.group_bin_files    = { "bacteria.map", "archaea.map" };
    cfg.reads              = { "virus.simulated.1.fq" };
    cfg.max_error          = 0;
    cfg.filter_hierarchy   = "1,2";
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_file ) );
}

SCENARIO( "Classify forced failure with different max. errors allowed", "[ganon-classify]" )
{
    auto cfg               = config_classify::defaultConfig();
    cfg.bloom_filter_files = { "bacteria.filter" };
    cfg.group_bin_files    = { "bacteria.map" };
    cfg.reads              = { "bacteria.simulated.1.fq" };
    cfg.max_error          = 0;

    const std::string undesired_output = "classify_output-b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_file, undesired_output ) );
}

SCENARIO( "Classify forced failure with different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg               = config_classify::defaultConfig();
    cfg.bloom_filter_files = { "bacteria.filter" };
    cfg.group_bin_files    = { "bacteria.map" };
    cfg.reads              = { "bacteria.simulated.1.fq" };
    cfg.max_error_unique   = 2;

    const std::string undesired_output = "classify_output-b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_file, undesired_output ) );
}


SCENARIO( "Classify multi-filter with partial matching reads", "[ganon-classify]" )
{

    // should match results as if working with only one filter because max_error = 0
    auto cfg1               = config_classify::defaultConfig();
    cfg1.bloom_filter_files = { "archaea.filter" };
    cfg1.group_bin_files    = { "archaea.map" };
    cfg1.reads              = { "archaea.simulated.1.fq" };
    cfg1.max_error          = 0;
    cfg1.output_file        = "a-a_e0.txt";

    auto cfg2               = config_classify::defaultConfig();
    cfg2.bloom_filter_files = { "bacteria.filter", "archaea.filter", "virus.filter" };
    cfg2.group_bin_files    = { "bacteria.map", "archaea.map", "virus.map" };
    cfg2.reads              = { "archaea.simulated.1.fq" };
    cfg2.max_error          = 0;
    cfg2.output_file        = "bav-a_e0.txt";

    REQUIRE( GanonClassify::run( cfg1 ) );
    REQUIRE( GanonClassify::run( cfg2 ) );
    REQUIRE( aux::filesAreEqual( cfg1.output_file, cfg2.output_file ) );
}

SCENARIO( "Classify multi-hierarchy with partial matching reads", "[ganon-classify]" )
{
    // should match results as if working with only one filter because max_error = 0
    // reads should 'survive' till the last filter
    auto cfg1               = config_classify::defaultConfig();
    cfg1.bloom_filter_files = { "virus.filter" };
    cfg1.group_bin_files    = { "virus.map" };
    cfg1.reads              = { "virus.simulated.1.fq" };
    cfg1.max_error          = 0;
    cfg1.output_file        = "v-v_e0.txt";


    // Test if the reads are surviving the hierachies
    auto cfg2               = config_classify::defaultConfig();
    cfg2.bloom_filter_files = { "bacteria.filter", "archaea.filter", "virus.filter" };
    cfg2.group_bin_files    = { "bacteria.map", "archaea.map", "virus.map" };
    cfg2.reads              = { "virus.simulated.1.fq" };
    cfg2.max_error          = 0;
    cfg2.filter_hierarchy   = "1,2,3";
    cfg2.output_file        = "bav-v_e0.txt";

    REQUIRE( GanonClassify::run( cfg1 ) );
    REQUIRE( GanonClassify::run( cfg2 ) );
    REQUIRE( aux::filesAreEqual( cfg1.output_file, cfg2.output_file ) );
}
