#include "CommandLineParser.hpp"

#include <defaults/defaults.hpp>

#include <cxxopts.hpp>

namespace GanonBuild
{

std::optional< Config > CommandLineParser::parse( int argc, char** argv )
{
    cxxopts::Options options( "ganon-build", "Ganon builder" );

    // clang-format off
    options.add_options()
        ( "i,input-file", "Define sequences to use. Tabular file with the fields: file [<tab> target]", cxxopts::value< std::string >() )
        ( "o,output-file", "Filter output file", cxxopts::value< std::string >() )
        ( "k,kmer-size", "k-mer size. Default: 19", cxxopts::value< uint8_t >() )
        ( "w,window-size", "window size. Default: 31", cxxopts::value< uint16_t >() )
        ( "s,hash-functions", "number of hash functions. 0 to auto-detect. Default: 0", cxxopts::value< uint8_t >() )
        ( "p,max-fp", "Maximum false positive rate per target. Used to define filter size [mutually exclusive --filter-size]. Default: 0.05", cxxopts::value< double >() )
        ( "f,filter-size", "Filter size (MB) [mutually exclusive --max-fp]", cxxopts::value< double >() )
        ( "j,mode", "mode to build filter [avg, smaller, smallest, faster, fastest]. Default: avg", cxxopts::value< std::string >() )
        ( "y,min-length", "min. sequence length (bp) to keep. 0 to keep all. Default: 0", cxxopts::value< uint64_t >() )
        ( "m,tmp-output-folder", "Folder to write temporary files", cxxopts::value< std::string >() )
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >())
        ( "verbose", "Verbose output mode", cxxopts::value<bool>())
        ( "quiet", "Quiet output mode", cxxopts::value<bool>())
        ( "h,help", "Show help commands" )
        ( "v,version", "Show current version" );
    // clang-format on

    const auto argcCopy = argc;
    const auto args     = options.parse( argc, argv );

    if ( argcCopy == 1 )
    {
        std::cerr << "Try 'ganon-build -h/--help' for more information." << std::endl;
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

    if ( args.count( "input-file" ) )
        config.input_file = args["input-file"].as< std::string >();
    if ( args.count( "output-file" ) )
        config.output_file = args["output-file"].as< std::string >();
    if ( args.count( "kmer-size" ) )
        config.kmer_size = args["kmer-size"].as< uint8_t >();
    if ( args.count( "window-size" ) )
        config.window_size = args["window-size"].as< uint16_t >();
    if ( args.count( "hash-functions" ) )
        config.hash_functions = args["hash-functions"].as< uint8_t >();
    if ( args.count( "max-fp" ) )
        config.max_fp = args["max-fp"].as< double >();
    if ( args.count( "filter-size" ) )
        config.filter_size = args["filter-size"].as< double >();
    if ( args.count( "mode" ) )
        config.mode = args["mode"].as< std::string >();
    if ( args.count( "min-length" ) )
        config.min_length = args["min-length"].as< uint64_t >();
    if ( args.count( "tmp-output-folder" ) )
        config.tmp_output_folder = args["tmp-output-folder"].as< std::string >();
    if ( args.count( "threads" ) )
        config.threads = args["threads"].as< uint16_t >();
    if ( args.count( "verbose" ) )
        config.verbose = args["verbose"].as< bool >();
    if ( args.count( "quiet" ) )
        config.quiet = args["quiet"].as< bool >();

    return config;
}

} // namespace GanonBuild
