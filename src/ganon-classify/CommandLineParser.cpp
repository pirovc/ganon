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
        ( "r,single-reads", "File[s] with single-end reads .fq .fastq .fasta .fa (e.g. file1.fq[.gz],[file2.fq[.gz] ... fileN.fq[.gz]])", cxxopts::value< std::vector< std::string > >() )
        ( "p,paired-reads", "Pairs of files with paired-end reads .fq .fastq .fasta .fa (e.g. file1.1.fq[.gz],file1.2.fq[.gz],[file2.1.fq[.gz],file2.2.fq[.gz] ... fileN.1.fq[.gz],fileN.2.fq[.gz]])", cxxopts::value< std::vector< std::string > >() )
        
        ( "i,ibf", "ibf (Interleaved Bloom Filter) file[s] (e.g. -b a.ibf,b.ibf OR -b a.ibf -b b.ibf )", cxxopts::value< std::vector< std::string > >() )
        ( "m,map", "map files[s]. Tab-separated file mapping target groups (taxids, assemblies) and bin identifiers with the following fields: target <tab> bin id (e.g. -g a.map,b.map OR -g a.map -g b.map)", cxxopts::value< std::vector< std::string > >() )
        ( "x,tax", "tax (taxonomy) files[s]. Tab-separated file with a complete tree with the following fields: node <tab> parent node <tab> rank <tab> name (e.g. -g a.tax,b.tax OR -g a.tax -g b.tax)", cxxopts::value< std::vector< std::string > >() )
        
        ( "c,hierarchy-labels", "Hierarchy labels for the database files (hierarchy follows the order of the sorted labels) (e.g. 1_host,2_target,1_host,3)", cxxopts::value< std::vector< std::string > >() )
        
        ( "k,min-kmers", "Minimum percentage of k-mers matching for a read to to be assigned [muttualy exclusive --max-error]. Default: 0.25", cxxopts::value< std::vector< float > >() )
        ( "e,max-error", "Maximum number of errors/mismatches allowed [muttualy exclusive --min-kmers]", cxxopts::value< std::vector< int16_t > >() )
        ( "u,max-error-unique", "Maximum number of errors/mismatches allowed for unique matches after filtering. Matches not passing this criterial will be flagged with negative k-mer counts.", cxxopts::value< std::vector< int16_t > >() )
        
        ( "f,offset", "Offset for skipping k-mers while counting. Function must be enabled on compilation time with -DGANON_OFFSET=ON. Default: 1 = no offset", cxxopts::value< uint16_t >() )
        
        ( "o,output-prefix", "Output prefix for output files (prefix.lca, prefix.rep, prefix.all, prefix.unclassified.fq). If multi-level hiearchy is provded, files are generated accordingly (prefix.hiearchy.lca, ...). Omit for output to STDOUT (only .lca will be printed)", cxxopts::value< std::string >() )
        ( "a,output-all", "Output file with all matches (prefix.all) [it can be very big]", cxxopts::value< bool >() )
        ( "n,output-unclassified-reads", "Output unclassified reads (prefix.unclassified.fq) file for unclassified reads", cxxopts::value< bool >() )
        ( "s,output-hierarchy-single", "Generate only one output (prefix.lca and prefix.rep) even with multiple hierarchy levels", cxxopts::value< bool >() )
        
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >())
        ( "n-batches", "Number of batches of n-reads to hold in memory. Default: 1000", cxxopts::value< uint32_t >())
        ( "n-reads", "Number of reads for each batch. Default: 400", cxxopts::value< uint32_t >())
        ( "verbose", "Verbose output mode", cxxopts::value< bool >())
        ( "quiet", "Quiet output mode (only outputs errors and warnings to the STDERR)", cxxopts::value< bool >())
        ( "h,help", "Print help" )
        ( "v,version", "Show version" );
    // clang-format on

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

    if ( args.count( "single-reads" ) )
        config.single_reads = args["single-reads"].as< std::vector< std::string > >();
    if ( args.count( "paired-reads" ) )
        config.paired_reads = args["paired-reads"].as< std::vector< std::string > >();

    if ( args.count( "ibf" ) )
        config.ibf = args["ibf"].as< std::vector< std::string > >();
    if ( args.count( "map" ) )
        config.map = args["map"].as< std::vector< std::string > >();
    if ( args.count( "tax" ) )
        config.tax = args["tax"].as< std::vector< std::string > >();
    if ( args.count( "hierarchy-labels" ) )
        config.hierarchy_labels = args["hierarchy-labels"].as< std::vector< std::string > >();

    if ( args.count( "min-kmers" ) )
        config.min_kmers = args["min-kmers"].as< std::vector< float > >();
    if ( args.count( "max-error" ) )
        config.max_error = args["max-error"].as< std::vector< int16_t > >();
    if ( args.count( "max-error-unique" ) )
        config.max_error_unique = args["max-error-unique"].as< std::vector< int16_t > >();
    if ( args.count( "offset" ) )
    {
#ifdef GANON_OFFSET
        config.offset = args["offset"].as< uint16_t >();
#else
        config.offset = 1;
#endif
    }

    if ( args.count( "output-prefix" ) )
        config.output_prefix = args["output-prefix"].as< std::string >();
    if ( args.count( "output-all" ) )
        config.output_all = args["output-all"].as< bool >();
    if ( args.count( "output-unclassified-reads" ) )
        config.output_unclassified_reads = args["output-unclassified-reads"].as< bool >();
    if ( args.count( "output-hierarchy-single" ) )
        config.output_hierarchy_single = args["output-hierarchy-single"].as< bool >();

    if ( args.count( "threads" ) )
        config.threads = args["threads"].as< uint16_t >();

    if ( args.count( "n-batches" ) )
        config.n_batches = args["n-batches"].as< uint32_t >();
    if ( args.count( "n-reads" ) )
        config.n_reads = args["n-reads"].as< uint32_t >();
    if ( args.count( "verbose" ) )
        config.verbose = args["verbose"].as< bool >();
    if ( args.count( "quiet" ) )
        config.quiet = args["quiet"].as< bool >();

    return config;
}

} // namespace GanonClassify
