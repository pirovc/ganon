#include <cxxopts.hpp>
#include <defaults/defaults.hpp>

struct Arguments
{

    static const uint64_t gbInBits = 8589934592;
    int                   argc;
    char**                argv;

    // Build options
    std::string seqid_bin_file;
    std::string output_file;
    uint64_t    bloom_filter_size;
    uint16_t    kmer_size;
    uint16_t    hash_functions;

    // Update options
    std::string update_bloom_filter_file;
    bool        update_complete;

    // General options
    std::vector< std::string > reference_files;
    uint16_t                   threads;
    bool verbose;


    Arguments( int _argc, char** _argv )
    : argc{ _argc }
    , argv{ _argv }
    {
    }

    bool parse()
    {

        cxxopts::Options options( "ganon-build", "Ganon builder" );

        // clang-format off
        options.add_options()
            ( "e,seqid-bin", "Tab-separated with the following fields: Seq. Identifier <tab> Pos. Seq. Start <tab> Pos. Seq. End <tab> Bin Id", cxxopts::value< std::string >() )
            ( "o,output-file", "Output file", cxxopts::value< std::string >() )
            ( "s,bloom-size", "Final bloom filter size in GB", cxxopts::value< uint64_t >()->default_value( "16" ) )
            ( "b,bloom-size-bits", "Final bloom filter size in bits", cxxopts::value< uint64_t >()->default_value( "0" ) )
            ( "k,kmer-size", "k size", cxxopts::value< uint16_t >()->default_value( "19" ) )
            ( "n,hash-functions", "Number of hash functions", cxxopts::value< uint16_t >()->default_value( "3" ) )
            ( "u,update-bloom-filter", "If provided, filte updated with new sequences", cxxopts::value< std::string >()->default_value( "" ) )
            ( "c,update-complete", "Old and new sequences are provided for updated bins", cxxopts::value< bool >()->default_value( "false" ) )
            ( "t,threads", "Number of threads", cxxopts::value< uint16_t >()->default_value( "1" ) )
            //( "silent", "Silent mode", cxxopts::value<bool>()->default_value("false"))
            ( "verbose", "Verbose output mode", cxxopts::value<bool>()->default_value("false"))
            ( "h,help", "Print help" )
            ( "v,version", "Show version" )
            ( "references", "references", cxxopts::value< std::vector< std::string > >() );
        // clang-format on
        options.parse_positional( { "references" } );
        options.positional_help( "reference.fna[.gz] [reference2.fna[.gz] ... referenceN.fna[.gz]]" );

        auto args = options.parse( argc, argv );

        if ( args.count( "version" ) )
        {
            std::cerr << "version: " << defaults::version_string << std::endl;
            return false;
        }
        else if ( args.count( "help" ) )
        {
            std::cerr << options.help() << std::endl;
            return false;
        }
        else
        {
            seqid_bin_file           = args["seqid-bin"].as< std::string >();
            update_bloom_filter_file = args["update-bloom-filter"].as< std::string >();
            update_complete          = args["update-complete"].as< bool >();
            threads                  = args["threads"].as< uint16_t >();
            output_file              = args["output-file"].as< std::string >();
            reference_files          = args["references"].as< std::vector< std::string > >();
            verbose          = args["verbose"].as< bool >();

            // Skip variables if updating
            if ( !update_bloom_filter_file.empty() )
            {
                std::cerr << "--bloom-size[-bits], --kmer-size --hash-funtions ignored, using metadata from "
                             "--update-bloom-filter file"
                          << std::endl;
            }
            else
            {
                kmer_size      = args["kmer-size"].as< uint16_t >();
                hash_functions = args["hash-functions"].as< uint16_t >();
                if ( args["bloom-size-bits"].as< uint64_t >() > 0 )
                    bloom_filter_size = args["bloom-size-bits"].as< uint64_t >();
                else
                    bloom_filter_size = args["bloom-size"].as< uint64_t >() * gbInBits;
            }
            return true;
        }
    }

    void print()
    {
        std::cerr << "seqid-bin: " << seqid_bin_file << std::endl;
        std::cerr << "bloom-size: " << std::fixed << std::setprecision( 2 )
                  << (float) bloom_filter_size / (float) gbInBits << std::endl;
        std::cerr << "bloom-size-bits: " << bloom_filter_size << std::endl;
        std::cerr << "kmer-size: " << kmer_size << std::endl;
        std::cerr << "hash-functions: " << hash_functions << std::endl;
        std::cerr << "bloom-filter: " << update_bloom_filter_file << std::endl;
        std::cerr << "update-complete: " << update_complete << std::endl;
        std::cerr << "threads: " << threads << std::endl;
        std::cerr << "verbose: " << verbose << std::endl;
        std::cerr << "output-file: " << output_file << std::endl;
        std::cerr << "references: " << std::endl;
        for ( const auto& reference_file : reference_files )
        {
            std::cerr << reference_file << std::endl;
        }
    }
};