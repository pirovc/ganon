#include "aux/Aux.hpp"

#include <ganon-classify/Config.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <catch2/catch.hpp>

namespace config_classify
{

GanonClassify::Config defaultConfig()
{
    GanonClassify::Config cfg;
    cfg.output_file      = "classify_output.txt";
    cfg.max_error        = 3;
    cfg.max_error_unique = -1;
    cfg.verbose          = true;
    cfg.testing          = true;
    cfg.threads          = 3;
    return cfg;
}

} // namespace config_classify

SCENARIO( "Classify", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter" };
    cfg.group_bin_files              = { "bacteria.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq" };
    const std::string desired_output = "b-b_e3.txt";

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

    const std::string desired_output   = "b-b_e3u1.txt";
    const std::string undesired_output = "b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_file, undesired_output ) );
}

SCENARIO( "Classify forced failure with different max. errors allowed", "[ganon-classify]" )
{
    auto cfg               = config_classify::defaultConfig();
    cfg.bloom_filter_files = { "bacteria.filter" };
    cfg.group_bin_files    = { "bacteria.map" };
    cfg.reads              = { "bacteria.simulated.1.fq" };
    cfg.max_error          = 0;

    const std::string undesired_output = "b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_file, undesired_output ) );
}


SCENARIO( "Classify multi-filter", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.bloom_filter_files           = { "bacteria.filter", "archaea.filter" };
    cfg.group_bin_files              = { "bacteria.map", "archaea.map" };
    cfg.reads                        = { "bacteria.simulated.1.fq", "archaea.simulated.1.fq" };
    const std::string desired_output = "ba-ba_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, desired_output ) );
}
