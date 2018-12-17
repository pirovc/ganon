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
        ( "c,filter-hierarchy", "Hierarchy of the given filter files (e.g. 1,1,2,3)", cxxopts::value< std::string >()->default_value( "" ) )
        ( "e,max-error", "Maximum number of errors/mismatches allowed for a match to be considered", cxxopts::value< int >()->default_value( "3" ) )
        ( "u,max-error-unique", "Maximum number of errors/mismatches allowed for unique matches after filtering. Matches not passing this criterial will have negative k-mer counts.", cxxopts::value< int >() )
        ( "o,output-file", "Output file with classification (omit for STDOUT). ", cxxopts::value< std::string >()->default_value( "" ) )
        ( "output-unclassified-file", "Output file for unclassified reads", cxxopts::value< std::string >()->default_value( "" ) )
        // option to skip filtering and work as a k-mer counter
        // option to output read len?
        ( "t,threads", "Number of threads", cxxopts::value< int >()->default_value( "3" ) )
        ( "verbose", "Verbose output mode", cxxopts::value<bool>()->default_value("false"))
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

    config.bloom_filter_files       = args["bloom-filter"].as< std::vector< std::string > >();
    config.group_bin_files          = args["group-bin"].as< std::vector< std::string > >();
    config.output_file              = args["output-file"].as< std::string >();
    config.output_unclassified_file = args["output-unclassified-file"].as< std::string >();
    config.max_error                = args["max-error"].as< int >();
    config.threads                  = args["threads"].as< int >();
    config.reads                    = args["reads"].as< std::vector< std::string > >();
    config.verbose                  = args["verbose"].as< bool >();
    config.filter_hierarchy         = args["filter-hierarchy"].as< std::string >();
    config.max_error_unique         = args.count( "max-error-unique" ) ? args["max-error-unique"].as< int >() : -1;

    return config;
}

} // namespace GanonClassify
