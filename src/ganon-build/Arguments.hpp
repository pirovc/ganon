#include <cxxopts.hpp>
#include <defaults/defaults.hpp>

struct Arguments
{

    static const uint64_t MBinBits = 8388608;
    // static const uint64_t GBinBits = 8589934592;

    int    argc;
    char** argv;

    // Build options
    std::string seqid_bin_file;
    std::string output_filter_file;
    uint64_t    filter_size;
    uint16_t    kmer_size;
    uint16_t    hash_functions;

    // Update options
    std::string update_filter_file;
    bool        update_complete;

    // General options
    std::vector< std::string > reference_files;
    uint16_t                   threads;
    uint16_t                   build_threads;
    bool                       verbose;


    Arguments( int _argc, char** _argv )
    : argc{ _argc }
    , argv{ _argv }
    {
    }

    bool parse()
    {

        cxxopts::Options options( "ganon-build", "Ganon builder" );

        int _argc = argc;

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
            seqid_bin_file     = args["seqid-bin-file"].as< std::string >();
            output_filter_file = args["output-filter-file"].as< std::string >();
            update_filter_file = args["update-filter-file"].as< std::string >();
            reference_files    = args["reference-files"].as< std::vector< std::string > >();
            update_complete    = args["update-complete"].as< bool >();
            threads            = args["threads"].as< uint16_t >();
            verbose            = args["verbose"].as< bool >();

            build_threads = threads - 1; //-1 reading files

            // Skip variables if updating, loads from existing filter file
            if ( !update_filter_file.empty() )
            {
                std::cerr << "--filter-size[-bits], --kmer-size --hash-funtions ignored, using metadata from "
                             "--update-filter-file"
                          << std::endl;
            }
            else
            {
                kmer_size      = args["kmer-size"].as< uint16_t >();
                hash_functions = args["hash-functions"].as< uint16_t >();
                if ( args["filter-size-bits"].as< uint64_t >() > 0 )
                    filter_size = args["filter-size-bits"].as< uint64_t >();
                else
                    filter_size = args["filter-size"].as< uint64_t >() * MBinBits;
            }
            return true;
        }
    }

    void print()
    {
        std::cerr << std::endl;
        std::cerr << "--seqid-bin-file      " << seqid_bin_file << std::endl;
        std::cerr << "--output-filter-file  " << output_filter_file << std::endl;
        std::cerr << "--update-filter-file  " << update_filter_file << std::endl;
        std::cerr << "--reference-files     " << std::endl;
        for ( const auto& reference_file : reference_files )
        {
            std::cerr << "                      " << reference_file << std::endl;
        }
        std::cerr << "--update-complete     " << update_complete << std::endl;
        std::cerr << "--filter-size         " << std::fixed << std::setprecision( 2 )
                  << (float) filter_size / (float) MBinBits << std::endl;
        std::cerr << "--filter-size-bits    " << filter_size << std::endl;
        std::cerr << "--kmer-size           " << kmer_size << std::endl;
        std::cerr << "--hash-functions      " << hash_functions << std::endl;
        std::cerr << "--threads             " << threads << std::endl;
        std::cerr << "--verbose             " << verbose << std::endl;
        std::cerr << std::endl;
    }
};