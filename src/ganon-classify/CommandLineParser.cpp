#include "CommandLineParser.hpp"

#include <defaults/defaults.hpp>

#include <cxxopts.hpp>

std::vector< std::string > split( const std::string& s, char delimiter )
{
    std::vector< std::string > tokens;
    std::string                token;
    std::istringstream         tokenStream( s );
    while ( std::getline( tokenStream, token, delimiter ) )
        tokens.push_back( token );
    return tokens;
}

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
        //( "silent", "Silent mode, just print results", cxxopts::value< int >()->default_value( "3" ) )
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

    config.clas_threads        = config.threads - 2; //-1 reading, -1 printing
    config.output_unclassified = false;
    config.unique_filtering    = false;
    config.max_error_unique    = 0;

    if ( !config.output_unclassified_file.empty() )
    {
        config.output_unclassified = true;
        config.clas_threads        = config.clas_threads - 1; //-1 printing unclassified
    }

    if ( args.count( "max-error-unique" ) )
    {
        config.max_error_unique = args["max-error-unique"].as< int >();
        if ( config.max_error_unique < config.max_error )
            config.unique_filtering = true;
    }


    if ( config.filter_hierarchy.empty() )
    {
        if ( config.bloom_filter_files.size() != config.group_bin_files.size() )
        {
            std::cerr << "Filter and group-bin files do not match" << std::endl;
            return std::nullopt;
        }
        else
        {
            for ( uint16_t h = 0; h < config.bloom_filter_files.size(); ++h )
            {
                config.filters["1"].push_back(
                    std::make_tuple( config.bloom_filter_files[h], config.group_bin_files[h] ) );
            }
        }
    }
    else
    {
        std::vector< std::string > hierarchy = split( config.filter_hierarchy, ',' );
        if ( hierarchy.size() != config.bloom_filter_files.size() || hierarchy.size() != config.group_bin_files.size() )
        {
            std::cerr << "Hierarchy does not match with the number of provided files" << std::endl;
            return std::nullopt;
        }
        else
        {
            for ( uint16_t h = 0; h < hierarchy.size(); ++h )
                config.filters[hierarchy[h]].push_back(
                    std::make_tuple( config.bloom_filter_files[h], config.group_bin_files[h] ) );
        }
    }
    return config;
}
