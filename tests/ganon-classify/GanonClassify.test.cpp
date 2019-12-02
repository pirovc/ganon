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
    cfg.output_prefix    = "classify_test_output";
    cfg.output_single    = true;
    cfg.threads          = 3;
    cfg.verbose          = false;
    cfg.quiet            = true;
    cfg.hierarchy_labels = { "1" };
    cfg.offset           = 1;

    return cfg;
}

} // namespace config_classify

// Static results

SCENARIO( "Classify", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter" };
    cfg.map                          = { "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error                    = { 3 };
    const std::string desired_output = "results/classify_output-b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify paired-reads concat", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter" };
    cfg.map                          = { "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax" };
    cfg.paired_reads                 = { "reads/bacteria_id.1.fq", "reads/bacteria_id.2.fq" };
    cfg.max_error                    = { 0 };
    const std::string desired_output = "results/classify_output-b-b_e0i1.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify paired-reads concat with unique errors", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter" };
    cfg.map                          = { "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax" };
    cfg.paired_reads                 = { "reads/bacteria_id.1.fq", "reads/bacteria_id.2.fq" };
    cfg.max_error                    = { 1 };
    cfg.max_error_unique             = { 0 };
    const std::string desired_output = "results/classify_output-b-b_e1u0i1.txt";

    REQUIRE( GanonClassify::run( cfg ) );

    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) == true );
}

SCENARIO( "Classify paired-reads and single-reads with multiple indices", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/archaea.filter", "filters/bacteria.filter" };
    cfg.map                          = { "files/archaea.map", "files/bacteria.map" };
    cfg.tax                          = { "files/archaea.tax", "files/bacteria.tax" };
    cfg.paired_reads                 = { "reads/bacteria_id.1.fq", "reads/bacteria_id.2.fq" };
    cfg.single_reads                 = { "reads/archaea.simulated.1.fq" };
    cfg.max_error                    = { 0 };
    cfg.hierarchy_labels             = { "1", "2" };
    const std::string desired_output = "results/classify_output-ab-ab_e0i1.txt";

    REQUIRE( GanonClassify::run( cfg ) );

    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) == true );
}

#ifdef GANON_OFFSET
SCENARIO( "Classify with offset", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter" };
    cfg.map                          = { "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq" };
    cfg.offset                       = 2;
    cfg.max_error                    = { 3 };
    const std::string desired_output = "results/classify_output-b-b_e3f2.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}
#endif

SCENARIO( "Classify with no errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter" };
    cfg.map                          = { "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error                    = { 0 };
    const std::string desired_output = "results/classify_output-b-b_e0.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify with min kmers", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter" };
    cfg.map                          = { "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq" };
    cfg.min_kmers                    = { 0.3 }; // should work the same as -e 3, threshold = 25 19-mers
    const std::string desired_output = "results/classify_output-b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify with different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter" };
    cfg.map                          = { "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error_unique             = { 1 };
    cfg.max_error                    = { 3 };
    const std::string desired_output = "results/classify_output-b-b_e3u1.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}


#ifdef GANON_OFFSET
SCENARIO( "Classify with offset and different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.filter" };
    cfg.map              = { "files/bacteria.map" };
    cfg.tax              = { "files/bacteria.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error        = { 2 };
    cfg.max_error_unique = { 0 };
    cfg.offset           = 6;

    const std::string desired_output = "results/classify_output-b-b_e2u0f6.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify with offset, min kmers and different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.filter" };
    cfg.map              = { "files/bacteria.map" };
    cfg.tax              = { "files/bacteria.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq" };
    cfg.min_kmers        = { 0.53 }; // ceil((100-19+1)*0.53) = 44 = ((100-19+1)-(2*19)) [2 errors allowed]
    cfg.max_error_unique = { 0 };
    cfg.offset           = 6;

    const std::string desired_output = "results/classify_output-b-b_e2u0f6.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}
#endif

SCENARIO( "Classify multi-filter without errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map                          = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error                    = { 0 };
    const std::string desired_output = "results/classify_output-ba-ba_e0.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify multi-filter with errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map                          = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error                    = { 4 };
    const std::string desired_output = "results/classify_output-ba-ba_e4.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify multi-filter with multiple errors", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map                          = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error                    = { 0, 4 };
    const std::string desired_output = "results/classify_output-ba-ba_e04.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify multi-hierarchy without errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map                          = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error                    = { 0 };
    cfg.hierarchy_labels             = { "1", "2" };
    const std::string desired_output = "results/classify_output-ba-ba_e0c12.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}

SCENARIO( "Classify multi-hierarchy split files", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map                          = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error                    = { 0 };
    cfg.hierarchy_labels             = { "1", "2" };
    std::string output_prefix1       = cfg.output_prefix + ".1.all";
    std::string output_prefix2       = cfg.output_prefix + ".2.all";
    cfg.output_single                = false;
    const std::string desired_output = "results/classify_output-ba-ba_e0c12.txt";

    REQUIRE( GanonClassify::run( cfg ) );

    int lines = aux::fileLines( output_prefix1 ) + aux::fileLines( output_prefix2 ) - 1;
    REQUIRE( lines == aux::fileLines( desired_output ) );
}

SCENARIO( "Classify multi-hierarchy with errors allowed", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map                          = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/bacteria.simulated.1.fq", "reads/archaea.simulated.1.fq" };
    cfg.max_error                    = { 4 };
    cfg.hierarchy_labels             = { "1", "2" };
    const std::string desired_output = "results/classify_output-ba-ba_e4c12.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}
SCENARIO( "Classify multi-hierarchy with multiple errors", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/archaea.filter", "filters/bacteria.filter" };
    cfg.map                          = { "files/archaea.map", "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/archaea.simulated.1.fq" };
    cfg.hierarchy_labels             = { "1", "2" };
    cfg.max_error                    = { 3, 4 };
    const std::string desired_output = "results/classify_output-ab-a_e34c12.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}
SCENARIO( "Classify multi-hierarchy with multiple errors and multiple unique errors", "[ganon-classify]" )
{
    auto cfg                         = config_classify::defaultConfig();
    cfg.ibf                          = { "filters/archaea.filter", "filters/bacteria.filter" };
    cfg.map                          = { "files/archaea.map", "files/bacteria.map" };
    cfg.tax                          = { "files/bacteria.tax", "files/archaea.tax" };
    cfg.single_reads                 = { "reads/archaea.simulated.1.fq" };
    cfg.hierarchy_labels             = { "1", "2" };
    cfg.max_error                    = { 3, 4 };
    cfg.max_error_unique             = { 0, 1 };
    const std::string desired_output = "results/classify_output-ab-a_e34c12u01.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::filesAreEqual( cfg.output_prefix + ".all", desired_output ) );
}


// Functionality

SCENARIO( "Classify problematic fastq", "[ganon-classify]" )
{
    auto cfg         = config_classify::defaultConfig();
    cfg.ibf          = { "filters/bacteria.filter" };
    cfg.map          = { "files/bacteria.map" };
    cfg.tax          = { "files/bacteria.tax" };
    cfg.single_reads = { "reads/problematic.fq" };
    cfg.max_error    = { 3 };

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileLines( cfg.output_prefix + ".all" ) == 5 );
}

SCENARIO( "Classify without matches", "[ganon-classify]" )
{
    auto cfg         = config_classify::defaultConfig();
    cfg.ibf          = { "filters/bacteria.filter" };
    cfg.map          = { "files/bacteria.map" };
    cfg.tax          = { "files/bacteria.tax" };
    cfg.single_reads = { "reads/virus.simulated.1.fq" };
    cfg.max_error    = { 0 };

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );
}

SCENARIO( "Classify multi-filter without matches", "[ganon-classify]" )
{
    auto cfg         = config_classify::defaultConfig();
    cfg.ibf          = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map          = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax          = { "files/bacteria.tax", "files/archaea.map" };
    cfg.single_reads = { "reads/virus.simulated.1.fq" };
    cfg.max_error    = { 0 };

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );
}

SCENARIO( "Classify multi-hierarchy without matches", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.filter", "filters/archaea.filter" };
    cfg.map              = { "files/bacteria.map", "files/archaea.map" };
    cfg.tax              = { "files/bacteria.tax", "files/archaea.map" };
    cfg.single_reads     = { "reads/virus.simulated.1.fq" };
    cfg.max_error        = { 0 };
    cfg.hierarchy_labels = { "1", "2" };

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );
}

SCENARIO( "Classify forced failure with different max. errors allowed", "[ganon-classify]" )
{
    auto cfg         = config_classify::defaultConfig();
    cfg.ibf          = { "filters/bacteria.filter" };
    cfg.map          = { "files/bacteria.map" };
    cfg.tax          = { "files/bacteria.tax" };
    cfg.single_reads = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error    = { 0 };

    const std::string undesired_output = "results/classify_output-b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_prefix + ".all", undesired_output ) );
}

SCENARIO( "Classify forced failure with different max. unique errors allowed", "[ganon-classify]" )
{
    auto cfg             = config_classify::defaultConfig();
    cfg.ibf              = { "filters/bacteria.filter" };
    cfg.map              = { "files/bacteria.map" };
    cfg.tax              = { "files/bacteria.tax" };
    cfg.single_reads     = { "reads/bacteria.simulated.1.fq" };
    cfg.max_error_unique = { 2 };
    cfg.max_error        = { 3 };

    const std::string undesired_output = "results/classify_output-b-b_e3.txt";

    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE_FALSE( aux::filesAreEqual( cfg.output_prefix + ".all", undesired_output ) );
}


SCENARIO( "Classify multi-filter with partial matching reads", "[ganon-classify]" )
{

    // should match results as if working with only one filter because max_error = 0
    auto cfg1          = config_classify::defaultConfig();
    cfg1.ibf           = { "filters/archaea.filter" };
    cfg1.map           = { "files/archaea.map" };
    cfg1.tax           = { "files/archaea.tax" };
    cfg1.single_reads  = { "reads/archaea.simulated.1.fq" };
    cfg1.max_error     = { 0 };
    cfg1.output_prefix = "a-a_e0";

    auto cfg2          = config_classify::defaultConfig();
    cfg2.ibf           = { "filters/bacteria.filter", "filters/archaea.filter", "filters/virus.filter" };
    cfg2.map           = { "files/bacteria.map", "files/archaea.map", "files/virus.map" };
    cfg2.tax           = { "files/bacteria.tax", "files/archaea.tax", "files/virus.tax" };
    cfg2.single_reads  = { "reads/archaea.simulated.1.fq" };
    cfg2.max_error     = { 0 };
    cfg2.output_prefix = "bav-a_e0";

    REQUIRE( GanonClassify::run( cfg1 ) );
    REQUIRE( GanonClassify::run( cfg2 ) );
    REQUIRE( aux::filesAreEqual( cfg1.output_prefix + ".all", cfg2.output_prefix + ".all" ) );
}

SCENARIO( "Classify multi-hierarchy with partial matching reads", "[ganon-classify]" )
{
    // should match results as if working with only one filter because max_error = 0
    // reads should 'survive' till the last filter
    auto cfg1          = config_classify::defaultConfig();
    cfg1.ibf           = { "filters/virus.filter" };
    cfg1.map           = { "files/virus.map" };
    cfg1.tax           = { "files/virus.tax" };
    cfg1.single_reads  = { "reads/virus.simulated.1.fq" };
    cfg1.max_error     = { 0 };
    cfg1.output_prefix = "v-v_e0";

    // Test if the reads are surviving the hierachies
    auto cfg2             = config_classify::defaultConfig();
    cfg2.ibf              = { "filters/bacteria.filter", "filters/archaea.filter", "filters/virus.filter" };
    cfg2.map              = { "files/bacteria.map", "files/archaea.map", "files/virus.map" };
    cfg2.tax              = { "files/bacteria.tax", "files/archaea.tax", "files/virus.tax" };
    cfg2.single_reads     = { "reads/virus.simulated.1.fq" };
    cfg2.max_error        = { 0 };
    cfg2.hierarchy_labels = { "1", "2", "3" };
    cfg2.output_prefix    = "bav-v_e0";

    REQUIRE( GanonClassify::run( cfg1 ) );
    REQUIRE( GanonClassify::run( cfg2 ) );
    REQUIRE( aux::filesAreEqual( cfg1.output_prefix + ".all", cfg2.output_prefix + ".all" ) );
}

SCENARIO( "Classify after update", "[ganon-classify]" )
{
    // No match
    auto cfg         = config_classify::defaultConfig();
    cfg.ibf          = { "filters/bacteria.filter" };
    cfg.map          = { "files/bacteria.map" };
    cfg.tax          = { "files/bacteria.tax" };
    cfg.single_reads = { "reads/virus.simulated.1.fq" };
    cfg.max_error    = { 3 };
    REQUIRE( GanonClassify::run( cfg ) );
    REQUIRE( aux::fileIsEmpty( cfg.output_prefix + ".all" ) );


    GanonBuild::Config cfg_build;
    cfg_build.update_filter_file = "filters/bacteria.filter";
    cfg_build.seqid_bin_file     = "files/bacteria_upd_virus_acc_bin.txt";
    cfg_build.output_filter_file = "bacteria_virus.filter";
    cfg_build.filter_size_bits   = 8388608;
    cfg_build.reference_files    = { "sequences/virus_NC_003676.1.fasta.gz",
                                  "sequences/virus_NC_011646.1.fasta.gz",
                                  "sequences/virus_NC_032412.1.fasta.gz",
                                  "sequences/virus_NC_035470.1.fasta.gz" };
    cfg_build.verbose            = cfg.verbose;
    cfg_build.quiet              = cfg.quiet;
    REQUIRE( GanonBuild::run( cfg_build ) );


    auto cfg2      = config_classify::defaultConfig();
    cfg2.max_error = { 3 };
    cfg2.ibf       = { cfg_build.output_filter_file };
    cfg2.map       = { "files/bacteria_virus.map" };
    cfg2.tax       = { "files/bacteria_virus.tax" };


    // Results of bacteria should be the same as without update (check if FP changed or filter was affected)
    cfg2.single_reads = { "reads/bacteria.simulated.1.fq" };
    REQUIRE( GanonClassify::run( cfg2 ) );
    REQUIRE( aux::filesAreEqual( cfg2.output_prefix + ".all", "results/classify_output-b-b_e3.txt" ) );

    // Virus reads should now map
    cfg2.single_reads = { "reads/virus.simulated.1.fq" };
    REQUIRE( GanonClassify::run( cfg2 ) );
    REQUIRE_FALSE( aux::fileIsEmpty( cfg2.output_prefix + ".all" ) );
}
