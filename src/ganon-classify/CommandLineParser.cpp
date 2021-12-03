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
        
        ( "i,ibf", "IBF (Interleaved Bloom Filter) file[s] generated with ganon-build (e.g. -b a.ibf,b.ibf OR -b a.ibf -b b.ibf )", cxxopts::value< std::vector< std::string > >() )
        ( "m,map", "Optional tab-separated file mapping bins ids (--ibf) to target groups/labels (e.g. taxids, assemblies). Targets can be repeated within/among filters if multiple bins represent the same group. If no --map file is provided, targets are bin ids. If multiple filters of hiearchies are provided, targets are: hierarchy label-filter id-bin id. Fields: target <tab> bin id (e.g. -g a.map,b.map OR -g a.map -g b.map)", cxxopts::value< std::vector< std::string > >() )
        ( "x,tax", "Optional tab-separated file with a taxonomic tree for LCA calculation. Will link targets provided in the --map files. Root node should be 1 with parent 0. Fields: node/target <tab> parent node <tab> rank <tab> name (e.g. -g a.tax,b.tax OR -g a.tax -g b.tax)", cxxopts::value< std::vector< std::string > >() )
        
        ( "y,hierarchy-labels", "Hierarchy labels to define level for database usage (hierarchy follows the order of the sorted labels) (e.g. 1_host,2_target,1_host,3). Default: 'H1'", cxxopts::value< std::vector< std::string > >() )
        
        ( "k,kmer-size", "k size to query - should be the same used to build filter. One per hierarchy label.", cxxopts::value< std::vector< uint8_t > >() )
        ( "w,window-size", "define window size for minimizers - should be the same used to build filter. One per hierarchy label.", cxxopts::value< std::vector< uint8_t > >() )
        ( "f,offset", "Offset for skipping k-mers while counting. One per hierarchy label. Default: 1 (no offset)", cxxopts::value< std::vector< uint8_t > >() )
        
        ( "c,rel-cutoff", "Relative cutoff (i.e. percentage of k-mers). 0 for no cutoff. One per filter [muttualy exclusive --rel-cutoff]. Default: 0.25", cxxopts::value< std::vector< double > >() )
        ( "b,abs-cutoff", "Absolute cutoff (i.e. number of errors). -1 for no cutoff. One per filter [muttualy exclusive --rel-cutoff].", cxxopts::value< std::vector< int16_t > >() )
        ( "d,rel-filter", "Relative filter. Additional percentage of matches allowed (relative to the best match). 1 for no filtering. One per hierarchy label [muttualy exclusive --abs-filter].", cxxopts::value< std::vector< double > >() )
        ( "e,abs-filter", "Absolute filter. Additional errors allowed (relative to the best match). -1 for no filtering. One per hierarchy label [muttualy exclusive --abs-filter]. Default: 0", cxxopts::value< std::vector< int16_t > >() )
        
        ( "o,output-prefix", "Output prefix (prefix.rep, [prefix.lca, prefix.all, prefix.unc]). If multi-level hierarchy is provided, files are generated accordingly (prefix.hierarchy.lca and prefix.hierarchy.all). Omit to output to STDOUT (only .rep will be printed)", cxxopts::value< std::string >() )
        ( "l,output-lca", "Output file with lca classification, one for each classified read (prefix.lca)", cxxopts::value< bool >() )
        ( "a,output-all", "Output file with all matches, one or more for each classified read (prefix.all) [it can be very big]", cxxopts::value< bool >() )
        ( "u,output-unclassified", "Output unclassified read ids (prefix.unc)", cxxopts::value< bool >() )
        ( "s,output-single", "Do not split output files (lca and all) with multiple hierarchy levels", cxxopts::value< bool >() )
        
        ( "t,threads", "Number of threads", cxxopts::value< uint16_t >())
        ( "n-reads", "Number of reads for each batch. Default: 400", cxxopts::value< uint32_t >())
        ( "n-batches", "Number of batches of n-reads to hold in memory. Default: 1000", cxxopts::value< uint32_t >())
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
    if ( args.count( "map" ) )
        config.map = args["map"].as< std::vector< std::string > >();
    if ( args.count( "tax" ) )
        config.tax = args["tax"].as< std::vector< std::string > >();
    if ( args.count( "hierarchy-labels" ) )
        config.hierarchy_labels = args["hierarchy-labels"].as< std::vector< std::string > >();

    if ( args.count( "kmer-size" ) )
        config.kmer_size = args["kmer-size"].as< std::vector< uint8_t > >();
    if ( args.count( "window-size" ) )
        config.window_size = args["window-size"].as< std::vector< uint8_t > >();
    if ( args.count( "offset" ) )
        config.offset = args["offset"].as< std::vector< uint8_t > >();

    if ( args.count( "rel-cutoff" ) )
        config.rel_cutoff = args["rel-cutoff"].as< std::vector< double > >();
    if ( args.count( "abs-cutoff" ) )
        config.abs_cutoff = args["abs-cutoff"].as< std::vector< int16_t > >();
    if ( args.count( "rel-filter" ) )
        config.rel_filter = args["rel-filter"].as< std::vector< double > >();
    if ( args.count( "abs-filter" ) )
        config.abs_filter = args["abs-filter"].as< std::vector< int16_t > >();

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

    if ( args.count( "threads" ) )
        config.threads = args["threads"].as< uint16_t >();
    if ( args.count( "n-reads" ) )
        config.n_reads = args["n-reads"].as< uint32_t >();
    if ( args.count( "n-batches" ) )
        config.n_batches = args["n-batches"].as< uint32_t >();
    if ( args.count( "verbose" ) )
        config.verbose = args["verbose"].as< bool >();
    if ( args.count( "quiet" ) )
        config.quiet = args["quiet"].as< bool >();

    return config;
}

} // namespace GanonClassify
