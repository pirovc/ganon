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
        ( "e,seqid-bin-file", "Tab-separated file linking sequences identifiers, start pos., end pos. and bin number. If not provided, iterate over input reference files and assumes one file=one bin. Sequence identifier is anything between the header start '>' and an empty space ' '. The file should contain the following fields: Seq. Identifier <tab> Pos. Seq. Start <tab> Pos. Seq. End <tab> Bin number", cxxopts::value< std::string >() )
        
        ( "o,output-filter-file", "Output file for filter (e.g. filter.ibf)", cxxopts::value< std::string >() )
        ( "u,update-filter-file", "Previously generated filter file to be updated", cxxopts::value< std::string >() )
        ( "c,update-complete", "When using --update-filter-file and all sequences are provided to update index, set this option to not only add sequences to the filter but also remove", cxxopts::value< bool >() )
        
        ( "f,false-positive", "False positive rate to build filter [mutually exclusive --bin-size-bits, --filter-size-mb]. Default: 0.05", cxxopts::value< double >() )
        ( "s,filter-size-mb", "Final filter size (MB) [mutually exclusive --bin-size-bits, --false-positive]", cxxopts::value< double >() )
        ( "b,bin-size-bits", "Bin size (bits) [mutually exclusive --filter-size-mb, --false-positive]", cxxopts::value< uint64_t >() )

        ( "k,kmer-size", "k-mer size to build filter (only forward strand). Default: 19", cxxopts::value< uint8_t >() )
        ( "w,window-size", "Window size. If set, filter is built with minimizers. ", cxxopts::value< uint32_t >() )
        ( "n,hash-functions", "Number of hash functions to build filter. Default: 3", cxxopts::value< uint16_t >() )
        ( "a,count-hashes", "Iterate over input to count the exact number of elements to insert into the filter", cxxopts::value<bool>())
        
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >())
        ( "n-refs", "Number of sequences for each batch. Default: 400", cxxopts::value< uint32_t >() )        
        ( "n-batches", "Number of batches of n-refs to hold in memory. Default: 1000", cxxopts::value< uint32_t >() )
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

    if ( args.count( "false-positive" ) )
        config.false_positive = args["false-positive"].as< double >();
    if ( args.count( "filter-size-mb" ) )
        config.filter_size_mb = args["filter-size-mb"].as< double >();
    if ( args.count( "bin-size-bits" ) )
        config.bin_size_bits = args["bin-size-bits"].as< uint64_t >();

    if ( args.count( "kmer-size" ) )
        config.kmer_size = args["kmer-size"].as< uint8_t >();
    if ( args.count( "window-size" ) )
        config.window_size = args["window-size"].as< uint32_t >();
    if ( args.count( "hash-functions" ) )
        config.hash_functions = args["hash-functions"].as< uint16_t >();
    if ( args.count( "count-hashes" ) )
        config.count_hashes = args["count-hashes"].as< bool >();

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
