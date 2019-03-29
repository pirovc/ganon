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
        ( "b,bloom-filter", "Input filter file[s]", cxxopts::value< std::vector< std::string > >() )
        ( "g,group-bin", "Tab-separated file[s] linking classification groups and bin identifiers. The file should contain the following fields: Group Identifier <tab> Bin Id", cxxopts::value< std::vector< std::string > >() )
        ( "c,filter-hierarchy", "Hierarchy of the given filter files (e.g. 1,1,2,3)", cxxopts::value< std::string >() )
        ( "e,max-error", "Maximum number of errors/mismatches allowed", cxxopts::value< std::string >() )
        ( "u,max-error-unique", "Maximum number of errors/mismatches allowed for unique matches after filtering. Matches not passing this criterial will be flagged with negative k-mer counts.", cxxopts::value< std::string >() )
        ( "m,min-kmers", "Minimum percentage of k-mers matching for a read to to be assigned [muttualy exclusive --max-error]", cxxopts::value< std::string >() )
        ( "f,offset", "Offset for skipping k-mers while counting. Function must be enabled on compilation time with -DGANON_OFFSET=ON. Default: 1 = no offset", cxxopts::value< uint16_t >() )
        ( "o,output-file", "Output file with classification (omit for STDOUT). ", cxxopts::value< std::string >() )
        ( "n,output-unclassified-file", "Output file for unclassified reads", cxxopts::value< std::string >() )
        ( "s,split-output-file-hierarchy", "Split output classification by hierarchy (filename will be outputfilename_hierachy", cxxopts::value< bool >() )
        ( "n-batches", "Number of batches of n-reads to hold in memory. Default: 1000", cxxopts::value< uint32_t >())
        ( "n-reads", "Number of reads for each batch. Default: 400", cxxopts::value< uint32_t >())
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >())
        ( "verbose", "Verbose output mode", cxxopts::value< bool >())
        ( "h,help", "Print help" )
        ( "v,version", "Show version" )
        ( "reads", "reads", cxxopts::value< std::vector< std::string > >() );
    // clang-format on


    options.parse_positional( { "reads" } );
    options.positional_help( "file1.fastq[.gz] [file2.fastq[.gz] ... fileN.fastq[.gz]]" );

    const auto argcCopy = argc;
    const auto args     = options.parse( argc, argv );

    if ( args.count( "help" ) || argcCopy == 1 )
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

    // Required
    config.bloom_filter_files = args["bloom-filter"].as< std::vector< std::string > >();
    config.group_bin_files    = args["group-bin"].as< std::vector< std::string > >();
    config.reads              = args["reads"].as< std::vector< std::string > >();

    // Default
    if ( args.count( "output-file" ) )
        config.output_file = args["output-file"].as< std::string >();
    if ( args.count( "output-unclassified-file" ) )
        config.output_unclassified_file = args["output-unclassified-file"].as< std::string >();
    if ( args.count( "max-error" ) )
        config.max_error = args["max-error"].as< std::string >();
    if ( args.count( "min-kmers" ) )
        config.min_kmers = args["min-kmers"].as< std::string >();
    if ( args.count( "offset" ) )
    {
#ifdef GANON_OFFSET
        config.offset = args["offset"].as< uint16_t >();
#else
        config.offset = 1;
#endif
    }
    if ( args.count( "threads" ) )
        config.threads = args["threads"].as< uint16_t >();
    if ( args.count( "split-output-file-hierarchy" ) )
        config.split_output_file_hierarchy = args["split-output-file-hierarchy"].as< bool >();
    if ( args.count( "verbose" ) )
        config.verbose = args["verbose"].as< bool >();
    if ( args.count( "filter-hierarchy" ) )
        config.filter_hierarchy = args["filter-hierarchy"].as< std::string >();
    if ( args.count( "max-error-unique" ) )
        config.max_error_unique = args["max-error-unique"].as< std::string >();
    if ( args.count( "n-batches" ) )
        config.n_batches = args["n-batches"].as< uint32_t >();
    if ( args.count( "n-reads" ) )
        config.n_reads = args["n-reads"].as< uint32_t >();

    return config;
}

} // namespace GanonClassify
