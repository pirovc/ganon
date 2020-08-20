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
        ( "r,reference-files", "Sequence files .fasta .fa .fna (e.g ref.fna[.gz],[ref2.fna[.gz],...,refN.fna[.gz]])", cxxopts::value< std::vector< std::string > >() )
        ( "d,directory-reference-files", "Directory with reference files. Do not provide wildcards, just path (e.g. /path/to/folder/)", cxxopts::value< std::string >() )
        ( "x,extension", "Extension of the files to search in the --directory-reference-files (e.g. '.fna')", cxxopts::value< std::string >() )

        ( "e,seqid-bin-file", "Tab-separated file linking sequences and bin identifiers. The file should contain the following fields: Seq. Identifier <tab> Pos. Seq. Start <tab> Pos. Seq. End <tab> Bin Id", cxxopts::value< std::string >() )
        ( "o,output-filter-file", "Output file for filter (e.g. filter.ibf)", cxxopts::value< std::string >() )
        ( "u,update-filter-file", "Previously generated filter file to be updated", cxxopts::value< std::string >() )
        ( "c,update-complete", "When using --update-filter-file and all sequences are provided to update index, set this option to not only add sequences to the filter but also remove", cxxopts::value< bool >() )
        
        ( "s,filter-size", "Final filter size in Megabytes (MB) [mutually exclusive --filter-size-bits]", cxxopts::value< uint32_t >() )
        ( "b,filter-size-bits", "Final filter size in Bits (bit) [mutually exclusive --filter-size]", cxxopts::value< uint64_t >() )
        
        ( "k,kmer-size", "k size to build filter", cxxopts::value< uint16_t >() )
        ( "n,hash-functions", "Number of hash functions to build filter", cxxopts::value< uint16_t >() )

        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >())
        ( "n-refs", "Number of sequences for each batch", cxxopts::value< uint32_t >() )        
        ( "n-batches", "Number of batches of n-refs to hold in memory", cxxopts::value< uint32_t >() )
        ( "verbose", "Verbose output mode", cxxopts::value<bool>())
        ( "quiet", "Quiet output mode (only outputs errors and warnings to the stderr)", cxxopts::value<bool>())
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

    if ( args.count( "reference-files" ) )
        config.reference_files = args["reference-files"].as< std::vector< std::string > >();
    if ( args.count( "directory-reference-files" ) )
        config.directory_reference_files = args["directory-reference-files"].as< std::string >();
    if ( args.count( "extension" ) )
        config.extension = args["extension"].as< std::string >();

    if ( args.count( "seqid-bin-file" ) )
        config.seqid_bin_file = args["seqid-bin-file"].as< std::string >();
    if ( args.count( "output-filter-file" ) )
        config.output_filter_file = args["output-filter-file"].as< std::string >();
    if ( args.count( "update-filter-file" ) )
        config.update_filter_file = args["update-filter-file"].as< std::string >();
    if ( args.count( "update-complete" ) )
        config.update_complete = args["update-complete"].as< bool >();

    if ( args.count( "filter-size" ) )
        config.filter_size = args["filter-size"].as< uint32_t >();
    if ( args.count( "filter-size-bits" ) )
        config.filter_size_bits = args["filter-size-bits"].as< uint64_t >();

    if ( args.count( "kmer-size" ) )
        config.kmer_size = args["kmer-size"].as< uint16_t >();
    if ( args.count( "hash-functions" ) )
        config.hash_functions = args["hash-functions"].as< uint16_t >();

    if ( args.count( "threads" ) )
        config.threads = args["threads"].as< uint16_t >();
    if ( args.count( "n-refs" ) )
        config.n_refs = args["n-refs"].as< uint32_t >();
    if ( args.count( "n-batches" ) )
        config.n_batches = args["n-batches"].as< uint32_t >();
    if ( args.count( "verbose" ) )
        config.verbose = args["verbose"].as< bool >();
    if ( args.count( "quiet" ) )
        config.quiet = args["quiet"].as< bool >();

    return config;
}

} // namespace GanonBuild
