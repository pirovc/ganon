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
        ( "e,seqid-bin-file", "Tab-separated file linking sequences and bin identifiers. The file should contain the following fields: Seq. Identifier <tab> Pos. Seq. Start <tab> Pos. Seq. End <tab> Bin Id", cxxopts::value< std::string >() )
        ( "o,output-filter-file", "Output filter file", cxxopts::value< std::string >() )
        ( "u,update-filter-file", "Previously generated filter file to be updated", cxxopts::value< std::string >()->default_value( "" ) )
        ( "c,update-complete", "Old and new sequences are provided for updated bins (used to remove sequences)", cxxopts::value< bool >()->default_value( "false" ) )
        ( "s,filter-size", "Final filter size in Megabytes (MB) [mutually exclusive --filter-size-bits]", cxxopts::value< uint64_t >()->default_value( "16" ) )
        ( "b,filter-size-bits", "Final filter size in Bits (bit) [mutually exclusive --filter-size]", cxxopts::value< uint64_t >()->default_value( "0" ) )
        ( "k,kmer-size", "k size", cxxopts::value< uint16_t >()->default_value( "19" ) )
        ( "n,hash-functions", "Number of hash functions", cxxopts::value< uint16_t >()->default_value( "3" ) )
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >()->default_value( "1" ) )
        ( "verbose", "Verbose output mode", cxxopts::value<bool>()->default_value("false"))
        ( "h,help", "Show help commands" )
        ( "v,version", "Show current version" )
        ( "reference-files", "reference-files", cxxopts::value< std::vector< std::string > >() );
    // clang-format on

    options.parse_positional( { "reference-files" } );
    options.positional_help( "ref.fna[.gz] [ref2.fna[.gz] ... refN.fna[.gz]]" );

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

    config.seqid_bin_file     = args["seqid-bin-file"].as< std::string >();
    config.output_filter_file = args["output-filter-file"].as< std::string >();
    config.update_filter_file = args["update-filter-file"].as< std::string >();
    config.reference_files    = args["reference-files"].as< std::vector< std::string > >();
    config.update_complete    = args["update-complete"].as< bool >();
    config.threads            = args["threads"].as< uint16_t >();
    config.verbose            = args["verbose"].as< bool >();

    config.build_threads = config.threads - 1; // -1 reading files

    // Skip variables if updating, loads from existing filter file
    if ( !config.update_filter_file.empty() )
    {
        std::cerr << "--filter-size[-bits], --kmer-size --hash-funtions ignored, using metadata from "
                     "--update-filter-file"
                  << std::endl;

        // TODO: what about all values set in the else block and not set here?
    }
    else
    {
        config.kmer_size      = args["kmer-size"].as< uint16_t >();
        config.hash_functions = args["hash-functions"].as< uint16_t >();

        config.filter_size = args["filter-size-bits"].as< uint64_t >() > 0
                                 ? args["filter-size-bits"].as< uint64_t >()
                                 : args["filter-size"].as< uint64_t >() * Config::MBinBits;
    }

    return config;
}

} // namespace GanonBuild
