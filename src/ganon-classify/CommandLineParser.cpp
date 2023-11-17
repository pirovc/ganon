#include "CommandLineParser.hpp"

#include <defaults/defaults.hpp>

#include <cxxopts.hpp>

namespace GanonClassify
{

std::optional< Config > CommandLineParser::parse( int argc, char** argv )
{
    cxxopts::Options options( "ganon-classify", "Ganon classifier" );

    // clang-format off
    options.add_options()
        ( "r,single-reads", "single-end reads file[s] (comma-separated, flat or gzipped)", cxxopts::value< std::vector< std::string > >() )
        ( "p,paired-reads", "paired-end reads file[s] (comma-separated, flat or gzipped)", cxxopts::value< std::vector< std::string > >() )
        
        ( "i,ibf", "ibf file[s] from ganon-build (comma-separated)", cxxopts::value< std::vector< std::string > >() )
        ( "x,tax", "tax file[s] from ganon-build for LCA calculation (comma-separated)", cxxopts::value< std::vector< std::string > >() )
        
        ( "y,hierarchy-labels", "Hierarchy labels to define level for classification. Hierarchy follows order of the sorted labels. Default: H1", cxxopts::value< std::vector< std::string > >() )
        
        ( "c,rel-cutoff", "Relative cutoff (i.e. percentage of minimizers). 0 for no cutoff. One or one per filter (comma-separated). Default: 0.2", cxxopts::value< std::vector< double > >() )
        ( "d,rel-filter", "Relative filter. Additional percentage of matches allowed (relative to the best match). 1 for no filtering. one or one per hierarchy label (comma-separated). Default: 0.0", cxxopts::value< std::vector< double > >() )
        ( "f,fpr-query", "Min. False positive for a query. 1 for no filtering. one or one per hierarchy label (comma-separated). Default: 1.0", cxxopts::value< std::vector< double > >() )
        
        ( "o,output-prefix", "Output prefix (prefix.rep, [prefix.one, prefix.all, prefix.unc]). If multi-level --hierarchy-labels is provided, files are generated accordingly (prefix.hierarchy.one and prefix.hierarchy.all). Omit to output to STDOUT (only .rep will be printed)", cxxopts::value< std::string >() )
        ( "l,output-lca", "Runs and outputs file with lca classification (prefix.one)", cxxopts::value< bool >() )
        ( "a,output-all", "Outputs file with all matches (prefix.all)", cxxopts::value< bool >() )
        ( "u,output-unclassified", "Outputs unclassified read ids (prefix.unc)", cxxopts::value< bool >() )
        ( "s,output-single", "Do not split output files (lca and all) with multi-level --hierarchy-labels", cxxopts::value< bool >() )
        
        ( "hibf", "Input is an Hierarchical IBF (.hibf) generated from raptor.", cxxopts::value< bool >())
        ( "skip-lca", "Skip LCA step.", cxxopts::value< bool >())
        ( "tax-root-node", "Define alternative root node for LCA. Default: 1", cxxopts::value< std::string >())
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >())
        ( "n-batches", "Number of batches of n-reads to hold in memory. Default: 1000", cxxopts::value< size_t >())
        ( "n-reads", "Number of reads for each batch. Default: 400", cxxopts::value< size_t >())
        ( "verbose", "Verbose output mode", cxxopts::value< bool >())
        ( "quiet", "Quiet output mode (only outputs errors and warnings to the STDERR)", cxxopts::value< bool >())
        ( "h,help", "Print help" )
        ( "v,version", "Show version" );
    // clang-format on

    const auto argcCopy = argc;
    const auto args     = options.parse( argc, argv );

    if ( argcCopy == 1 )
    {
        std::cerr << "Try 'ganon-classify -h/--help' for more information." << std::endl;
        return std::nullopt;
    }
    else if ( args.count( "help" ) )
    {
        std::cerr << options.help() << std::endl;
        return std::nullopt;
    }
    else if ( args.count( "version" ) )
    {
        std::cerr << "version: " << defaults::version_string << std::endl;
        return std::nullopt;
    }

    Config config;

    if ( args.count( "single-reads" ) )
        config.single_reads = args["single-reads"].as< std::vector< std::string > >();
    if ( args.count( "paired-reads" ) )
        config.paired_reads = args["paired-reads"].as< std::vector< std::string > >();

    if ( args.count( "ibf" ) )
        config.ibf = args["ibf"].as< std::vector< std::string > >();
    if ( args.count( "tax" ) )
        config.tax = args["tax"].as< std::vector< std::string > >();

    if ( args.count( "hierarchy-labels" ) )
        config.hierarchy_labels = args["hierarchy-labels"].as< std::vector< std::string > >();

    if ( args.count( "rel-cutoff" ) )
        config.rel_cutoff = args["rel-cutoff"].as< std::vector< double > >();
    if ( args.count( "rel-filter" ) )
        config.rel_filter = args["rel-filter"].as< std::vector< double > >();
    if ( args.count( "fpr-query" ) )
        config.fpr_query = args["fpr-query"].as< std::vector< double > >();

    if ( args.count( "output-prefix" ) )
        config.output_prefix = args["output-prefix"].as< std::string >();
    if ( args.count( "output-lca" ) )
        config.output_lca = args["output-lca"].as< bool >();
    if ( args.count( "output-all" ) )
        config.output_all = args["output-all"].as< bool >();
    if ( args.count( "output-unclassified" ) )
        config.output_unclassified = args["output-unclassified"].as< bool >();
    if ( args.count( "output-single" ) )
        config.output_single = args["output-single"].as< bool >();

    if ( args.count( "hibf" ) )
        config.hibf = args["hibf"].as< bool >();
    if ( args.count( "skip-lca" ) )
        config.skip_lca = args["skip-lca"].as< bool >();
    if ( args.count( "tax-root-node" ) )
        config.tax_root_node = args["tax-root-node"].as< std::string >();
    if ( args.count( "threads" ) )
        config.threads = args["threads"].as< uint16_t >();
    if ( args.count( "n-reads" ) )
        config.n_reads = args["n-reads"].as< size_t >();
    if ( args.count( "n-batches" ) )
        config.n_batches = args["n-batches"].as< size_t >();
    if ( args.count( "verbose" ) )
        config.verbose = args["verbose"].as< bool >();
    if ( args.count( "quiet" ) )
        config.quiet = args["quiet"].as< bool >();

    return config;
}

} // namespace GanonClassify
