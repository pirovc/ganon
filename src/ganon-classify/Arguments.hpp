#include <cxxopts.hpp>
#include <defaults/defaults.hpp>

std::vector< std::string > split( const std::string& s, char delimiter )
{
    std::vector< std::string > tokens;
    std::string                token;
    std::istringstream         tokenStream( s );
    while ( std::getline( tokenStream, token, delimiter ) )
        tokens.push_back( token );
    return tokens;
}

struct Arguments
{

    int    argc;
    char** argv;

    std::string                                                                    output_file;
    std::string                                                                    output_unclassified_file;
    std::string                                                                    filter_hierarchy;
    uint16_t                                                                       max_error;
    uint16_t                                                                       threads;
    uint16_t                                                                       clas_threads;
    bool                                                                           output_unclassified;
    bool                                                                           unique_filtering;
    uint16_t                                                                       max_error_unique;
    std::vector< std::string >                                                     bloom_filter_files;
    std::vector< std::string >                                                     group_bin_files;
    std::vector< std::string >                                                     reads;
    std::map< std::string, std::vector< std::tuple< std::string, std::string > > > filters;
    bool                                                                           verbose;

    Arguments( int _argc, char** _argv )
    : argc{ _argc }
    , argv{ _argv }
    {
    }

    bool parse()
    {

        cxxopts::Options options( "ganon-classify", "Ganon classifier" );

        int _argc = argc;

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

        auto args = options.parse( argc, argv );

        if ( args.count( "help" ) || _argc == 1 )
        {
            std::cerr << options.help() << std::endl;
            return false;
        }
        else if ( args.count( "version" ) )
        {
            std::cerr << "version: " << defaults::version_string << std::endl;
            return false;
        }
        else
        {

            bloom_filter_files       = args["bloom-filter"].as< std::vector< std::string > >();
            group_bin_files          = args["group-bin"].as< std::vector< std::string > >();
            output_file              = args["output-file"].as< std::string >();
            output_unclassified_file = args["output-unclassified-file"].as< std::string >();
            max_error                = args["max-error"].as< int >();
            threads                  = args["threads"].as< int >();
            clas_threads             = threads - 2; //-1 reading, -1 printing
            output_unclassified      = false;
            unique_filtering         = false;
            reads                    = args["reads"].as< std::vector< std::string > >();
            verbose                  = args["verbose"].as< bool >();
            filter_hierarchy         = args["filter-hierarchy"].as< std::string >();

            if ( !output_unclassified_file.empty() )
            {
                output_unclassified = true;
                clas_threads        = clas_threads - 1; //-1 printing unclassified
            }

            if ( args.count( "max-error-unique" ) )
            {
                max_error_unique = args["max-error-unique"].as< int >();
                if ( max_error_unique < max_error )
                    unique_filtering = true;
            }
            else
            {
                max_error_unique = 0;
            }


            if ( filter_hierarchy.empty() )
            {
                if ( bloom_filter_files.size() != group_bin_files.size() )
                {
                    std::cerr << "Filter and group-bin files do not match" << std::endl;
                    return false;
                }
                else
                {
                    for ( uint16_t h = 0; h < bloom_filter_files.size(); ++h )
                    {
                        filters["1"].push_back( std::make_tuple( bloom_filter_files[h], group_bin_files[h] ) );
                    }
                }
            }
            else
            {
                std::vector< std::string > hierarchy = split( filter_hierarchy, ',' );
                if ( hierarchy.size() != bloom_filter_files.size() || hierarchy.size() != group_bin_files.size() )
                {
                    std::cerr << "Hierarchy does not match with the number of provided files" << std::endl;
                    return false;
                }
                else
                {
                    for ( uint16_t h = 0; h < hierarchy.size(); ++h )
                        filters[hierarchy[h]].push_back( std::make_tuple( bloom_filter_files[h], group_bin_files[h] ) );
                }
            }
            return true;
        }
    }

    void print()
    {
        std::cerr << std::endl;
        std::cerr << "--filter-hierarchy          " << filter_hierarchy << std::endl;
        std::cerr << "--bloom-filter,--group-bin  " << std::endl;
        for ( auto const& hierarchy : filters )
        {
            for ( auto const& file : hierarchy.second )
            {
                std::cerr << "                            (" << hierarchy.first << ") " << std::get< 0 >( file ) << ", "
                          << std::get< 1 >( file ) << std::endl;
            }
        }
        std::cerr << "--max-error                 " << max_error << std::endl;
        std::cerr << "--max-error-unique          " << max_error_unique << std::endl;
        std::cerr << "--output-file               " << output_file << std::endl;
        std::cerr << "--output-unclassified-file  " << output_unclassified_file << std::endl;
        std::cerr << "--verbose                   " << verbose << std::endl;
        std::cerr << "--threads                   " << threads << std::endl;
        std::cerr << "--reads                     " << std::endl;
        for ( const auto& s : reads )
            std::cerr << "                            " << s << std::endl;
        std::cerr << std::endl;
    }
};