#include "aux/Aux.hpp"

#include <ganon-classify/Config.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <catch2/catch.hpp>

namespace config_classify
{

Config testConfig()
{
    Config cfg;
    cfg.bloom_filter_files = { "classify.filter" };
    cfg.group_bin_files    = { "bacteria_group_bin.txt" };
    cfg.output_file        = "test_output.txt";
    cfg.max_error          = 3;
    cfg.reads              = { "simulated.1.fq" };
    cfg.verbose            = false;
    cfg.testing            = true;
    return cfg;
}

const std::string outputFile = "classify_output.txt";

} // namespace config_classify

SCENARIO( "Classify", "[ganon-classify]" )
{
    const auto cfg = config_classify::testConfig();
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_file, config_classify::outputFile ) );
}

SCENARIO( "Classify forced failure with different number of errors allowed", "[ganon-classify]" )
{
    auto cfg      = config_classify::testConfig();
    cfg.max_error = 1;

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_file, config_classify::outputFile ) );
}
