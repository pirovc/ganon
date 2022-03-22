#include "aux/Aux.hpp"

#include <seqan3/core/debug_stream.hpp>

#include <ganon-classify/Config.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

#include <filesystem>

using namespace seqan3::literals;

namespace config_classify
{

GanonClassify::Config defaultConfig( const std::string prefix )
{
    GanonClassify::Config cfg;
    cfg.output_prefix       = prefix;
    cfg.output_all          = true;
    cfg.output_lca          = true;
    cfg.output_unclassified = true;
    cfg.threads             = 4;
    cfg.kmer_size           = { 10 };
    cfg.verbose             = false;
    cfg.quiet               = true;
    return cfg;
}


typedef std::map< std::string, std::map< std::string, uint32_t > > TOut;
typedef std::vector< std::string >                                 TUnc;
typedef std::map< std::string, std::string >                       TMap;
typedef std::map< std::string, std::string >                       TTax;

struct Res
{

    void parse_rep( std::string file )
    {
        std::string   line;
        std::ifstream infile{ file };
        while ( std::getline( infile, line, '\n' ) )
        {
            std::istringstream         stream_line( line );
            std::vector< std::string > fields;
            std::string                field;
            while ( std::getline( stream_line, field, '\t' ) )
                fields.push_back( field );
            /* 0 hierarchy_label <tab>
               1 target <tab>
               2 matches <tab>
               3 unique reads <tab>
               4 lca reads <tab>
               5 [rank <tab>]
               6 [name]
            */
            if ( fields[0] == "#total_classified" )
            {
                total_classified = std::stoul( fields[1] );
            }
            else if ( fields[0] == "#total_unclassified" )
            {
                total_unclassified = std::stoul( fields[1] );
            }
            else
            {
                matches += std::stoul( fields[2] );
                unique_reads += std::stoul( fields[3] );
                lca_reads += std::stoul( fields[4] );
            }
        }
        total_reads = total_classified + total_unclassified;
        infile.close();
    }

    TOut parse_all_lca( std::string file, uint64_t& nlines )
    {
        TOut          out;
        std::string   line;
        std::ifstream infile{ file };
        while ( std::getline( infile, line, '\n' ) )
        {
            std::istringstream         stream_line( line );
            std::vector< std::string > fields;
            std::string                field;
            while ( std::getline( stream_line, field, '\t' ) )
                fields.push_back( field );
            // readid <tab> target <tab> k-mer count
            out[fields[0]][fields[1]] = std::stoul( fields[2] );
            nlines++;
        }
        infile.close();
        return out;
    }

    TUnc parse_unc( std::string file )
    {
        TUnc          out;
        std::string   line;
        std::ifstream infile{ file };
        // Read the next line from File untill it reaches the end.
        while ( std::getline( infile, line ) )
        {
            out.push_back( line );
        }
        return out;
    }

    Res()
    {
    }

    Res( GanonClassify::Config& cfg )
    {
        parse_rep( cfg.output_prefix + ".rep" );
        if ( cfg.output_all )
        {
            all = parse_all_lca( cfg.output_prefix + ".all", lines_all );
        }
        if ( cfg.output_lca )
        {
            lca = parse_all_lca( cfg.output_prefix + ".lca", lines_lca );
        }
        if ( cfg.output_unclassified )
        {
            unc = parse_unc( cfg.output_prefix + ".unc" );
        }
    }

    TOut     all;
    TOut     lca;
    TUnc     unc;
    uint64_t total_classified   = 0;
    uint64_t total_unclassified = 0;
    uint64_t total_reads        = 0;
    uint64_t matches            = 0;
    uint64_t unique_reads       = 0;
    uint64_t lca_reads          = 0;
    uint64_t lines_all          = 0;
    uint64_t lines_lca          = 0;
};

void sanity_check( GanonClassify::Config const& cfg, Res const& res )
{
    if ( cfg.output_all )
    {
        // Output as many reads and matches as reported (.rep)
        REQUIRE( res.all.size() == res.total_classified );
        REQUIRE( res.lines_all == res.matches );
    }

    if ( cfg.output_lca && !cfg.tax.empty() )
    {
        // Output as many lca reads as reported (.rep)
        REQUIRE( res.lca.size() == res.total_classified );
        REQUIRE( res.lines_lca == res.total_classified );
    }

    if ( cfg.output_unclassified )
    {
        // Output as many unclassified reads as reported (.rep)
        REQUIRE( res.unc.size() == res.total_unclassified );
    }
}

void write_map( std::string file, TMap const& map_data )
{
    std::ofstream outfile( file );
    for ( auto const& [binid, target] : map_data )
    {
        outfile << target << '\t' << binid << '\n';
    }
    outfile.close();
}

void write_tax( std::string file, TTax const& tax_data )
{
    // Automatically create name and ranks
    std::ofstream outfile( file );
    // Add root "1"
    outfile << "1" << '\t' << "0" << '\t' << "root" << '\t' << "root" << '\n';
    for ( auto const& [target, parent] : tax_data )
    {
        outfile << target << '\t' << parent << '\t' << "rank-" + target << '\t' << "name-" + target << '\n';
    }
    outfile.close();
}

} // namespace config_classify


SCENARIO( "classifying reads without errors", "[ganon-classify][without-errors]" )
{

    std::string folder_prefix{ "ganon-classify-wo-errors/" };
    std::filesystem::create_directory( folder_prefix );

    // Reads (14bp)
    aux::write_sequences( folder_prefix + "rA.fasta", { "AAAAAAAAAAAAAA"_dna4 }, { "readA" } );
    aux::write_sequences( folder_prefix + "rC.fasta", { "CCCCCCCCCCCCCC"_dna4 }, { "readC" } );
    aux::write_sequences( folder_prefix + "rT.fasta", { "TTTTTTTTTTTTTT"_dna4 }, { "readT" } );
    aux::write_sequences( folder_prefix + "rG.fasta", { "GGGGGGGGGGGGGG"_dna4 }, { "readG" } );

    // Sequences (20bp)
    const ids_type       ids{ "seqA", "seqC", "seqT", "seqG" };
    const sequences_type seqs{ "AAAAAAAAAAAAAAAAAAAA"_dna4,
                               "CCCCCCCCCCCCCCCCCCCC"_dna4,
                               "TTTTTTTTTTTTTTTTTTTT"_dna4,
                               "GGGGGGGGGGGGGGGGGGGG"_dna4 };


    std::string        base_prefix{ folder_prefix + "base_build1" };
    GanonBuild::Config cfg_build;
    cfg_build.bin_size_bits      = 5000;
    cfg_build.verbose            = false;
    cfg_build.quiet              = true;
    cfg_build.kmer_size          = 10;
    cfg_build.output_filter_file = base_prefix + ".ibf";
    cfg_build.reference_files    = aux::write_sequences_files( base_prefix, "fasta", seqs, ids );
    REQUIRE( GanonBuild::run( cfg_build ) );

    SECTION( "default params" )
    {
        std::string prefix{ folder_prefix + "default" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "rA.fasta" };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches on binids (Rev. Comp.)
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["0"] == 5 );
        REQUIRE( res.all["readA"]["2"] == 5 );

        SECTION( "without --output-lca" )
        {
            std::string prefix{ folder_prefix + "default_wo_output_lca" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.output_lca   = false;
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.single_reads = { folder_prefix + "rA.fasta" };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );
            REQUIRE_FALSE( std::filesystem::exists( prefix + ".lca" ) );
        }
        SECTION( "without --output-all" )
        {
            std::string prefix{ folder_prefix + "default_wo_output_all" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.output_all   = false;
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.single_reads = { folder_prefix + "rA.fasta" };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );
            REQUIRE_FALSE( std::filesystem::exists( prefix + ".all" ) );
        }
    }

    SECTION( "--paired-reads" )
    {
        std::string prefix{ folder_prefix + "paired_reads" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.paired_reads = { folder_prefix + "rA.fasta", folder_prefix + "rT.fasta" };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // will report header of first pair "readA" and match the rev.comp. of readT (and opposite)
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["0"] == 10 );
        REQUIRE( res.all["readA"]["2"] == 10 );
    }

    SECTION( "--single-reads and --paired-reads" )
    {
        std::string prefix{ folder_prefix + "single_paired_reads" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "rC.fasta", folder_prefix + "rG.fasta" };
        cfg.paired_reads = { folder_prefix + "rA.fasta", folder_prefix + "rT.fasta" };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // same as paired and single C and G
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readC"].size() == 2 );
        REQUIRE( res.all["readG"].size() == 2 );
        REQUIRE( res.all["readA"]["0"] == 10 );
        REQUIRE( res.all["readC"]["1"] == 5 );
        REQUIRE( res.all["readG"]["1"] == 5 );
        REQUIRE( res.all["readA"]["2"] == 10 );
        REQUIRE( res.all["readC"]["3"] == 5 );
        REQUIRE( res.all["readG"]["3"] == 5 );
    }

    SECTION( "--map" )
    {
        std::string prefix{ folder_prefix + "map" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "rA.fasta" };
        config_classify::write_map( prefix + ".map",
                                    { { "0", "AorT" }, { "1", "CorG" }, { "2", "AorT" }, { "3", "CorG" } } );
        cfg.map = { prefix + ".map" };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches on target (max. k-mer count)
        REQUIRE( res.all["readA"].size() == 1 );
        REQUIRE( res.all["readA"]["AorT"] == 5 );
    }

    SECTION( "--map and --tax" )
    {
        std::string prefix{ folder_prefix + "map_tax" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "rA.fasta" };
        config_classify::write_map( prefix + ".map", { { "0", "A" }, { "1", "C" }, { "2", "T" }, { "3", "G" } } );
        cfg.map = { prefix + ".map" };

        //     1
        //    ATCG
        //  AT    CG
        // A  T  C  G
        config_classify::write_tax( prefix + ".tax",
                                    { { "A", "AT" },
                                      { "C", "CG" },
                                      { "T", "AT" },
                                      { "G", "CG" },
                                      { "CG", "ATCG" },
                                      { "AT", "ATCG" },
                                      { "ATCG", "1" } } );
        cfg.tax = { prefix + ".tax" };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // All matches on targets from map
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["A"] == 5 );
        REQUIRE( res.all["readA"]["T"] == 5 );

        // LCA matches from tax
        REQUIRE( res.lca["readA"].size() == 1 );
        REQUIRE( res.lca["readA"]["AT"] == 5 );
    }

    SECTION( "incomplete --map and --tax" )
    {
        std::string prefix{ folder_prefix + "incomplete_map_tax" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "rA.fasta" };
        config_classify::write_map( prefix + ".map",
                                    { { "1", "C" }, { "2", "T" }, { "3", "G" } } ); // missing { "0", "A" }
        cfg.map = { prefix + ".map" };

        //     1
        //    ATCG
        //  AT    CG
        // A  T    G
        config_classify::write_tax( prefix + ".tax",
                                    { { "A", "AT" },
                                      { "T", "AT" },
                                      { "G", "CG" },
                                      { "CG", "ATCG" },
                                      { "AT", "ATCG" },
                                      { "ATCG", "1" } } ); // missing  { "C", "CG" }
        cfg.tax = { prefix + ".tax" };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // All matches on targets from map
        REQUIRE( res.all["readA"].size() == 1 );
        REQUIRE( res.all["readA"]["A"] == 0 ); // A not found
        REQUIRE( res.all["readA"]["T"] == 5 );

        // LCA matches  are the same as all (no A matchces)
        REQUIRE( res.lca["readA"].size() == 1 );
        REQUIRE( res.lca["readA"]["T"] == 5 );
    }

    SECTION( "--window-size" )
    {
        // build with --window-size
        std::string        base_prefix_ws{ folder_prefix + "base_build_window_size" };
        GanonBuild::Config cfg_build_ws;
        cfg_build_ws.bin_size_bits      = 5000;
        cfg_build_ws.quiet              = true;
        cfg_build_ws.kmer_size          = 10;
        cfg_build_ws.window_size        = 12;
        cfg_build_ws.output_filter_file = base_prefix_ws + ".ibf";
        cfg_build_ws.reference_files    = aux::write_sequences_files( base_prefix_ws, "fasta", seqs, ids );
        REQUIRE( GanonBuild::run( cfg_build_ws ) );

        std::string prefix{ folder_prefix + "window_size" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix_ws + ".ibf" };
        cfg.single_reads = { folder_prefix + "rA.fasta" };
        cfg.kmer_size    = { 10 };
        cfg.window_size  = { 12 };
        cfg.rel_filter   = { 1 };
        cfg.rel_cutoff   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // one kmer counted since they share same hash
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["0"] == 1 );
        REQUIRE( res.all["readA"]["2"] == 1 );

        SECTION( "--paired-reads" )
        {
            std::string prefix{ folder_prefix + "window_size_paired_reads" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix_ws + ".ibf" };
            cfg.paired_reads = { folder_prefix + "rA.fasta", folder_prefix + "rT.fasta" };
            cfg.kmer_size    = { 10 };
            cfg.window_size  = { 12 };
            cfg.rel_filter   = { 1 };
            cfg.rel_cutoff   = { 0 };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // will report header of first pair "readA" and match the rev.comp. of readT (and opposite)
            REQUIRE( res.all["readA"].size() == 2 );
            REQUIRE( res.all["readA"]["0"] == 2 );
            REQUIRE( res.all["readA"]["2"] == 2 );
        }
    }

    SECTION( "--offset 3" )
    {
        std::string prefix{ folder_prefix + "offset" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "rA.fasta" };
        cfg.offset       = { 3 };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches on binids (Rev. Comp.)
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["0"] == 3 );
        REQUIRE( res.all["readA"]["2"] == 3 );

        SECTION( "--paired-reads" )
        {
            std::string prefix{ folder_prefix + "offset_paired_reads" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.paired_reads = { folder_prefix + "rA.fasta", folder_prefix + "rT.fasta" };
            cfg.offset       = { 3 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Matches on binids (Rev. Comp.)
            REQUIRE( res.all["readA"].size() == 2 );
            REQUIRE( res.all["readA"]["0"] == 6 );
            REQUIRE( res.all["readA"]["2"] == 6 );
        }
    }

    SECTION( "wrong --kmer-size" )
    {
        std::string prefix{ folder_prefix + "wrong_kmer_size" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "rA.fasta" };
        cfg.kmer_size    = { 19 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );
        // do not classify and return empty file
        REQUIRE( res.all.empty() );
        // report unclassified
        REQUIRE( res.unc.size() == 1 );
        REQUIRE( res.unc[0] == "readA" );
    }

    SECTION( "without --output-prefix" )
    {
        // avoid printing on test results (ignoring STDOUT)
        std::streambuf* old = std::cout.rdbuf( 0 );

        std::string prefix{ folder_prefix + "wo_output_prefix" };
        auto        cfg   = config_classify::defaultConfig( prefix );
        cfg.ibf           = { base_prefix + ".ibf" };
        cfg.single_reads  = { folder_prefix + "rA.fasta" };
        cfg.output_prefix = "";
        REQUIRE( GanonClassify::run( cfg ) );
        // No files created
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.rep" ) );
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.lca" ) );
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.all" ) );
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.unc" ) );

        std::cout.rdbuf( old );
    }

    SECTION( "--hierarchy-labels" )
    {
        // Additional filter with repeated sequences (As) and new sequence (CG)

        // Reads (14bp)
        aux::write_sequences( folder_prefix + "rCG.fasta", { "CGCGCGCGCGCGCG"_dna4 }, { "readCG" } );

        // Sequences (20bp)
        const ids_type       ids2{ "seqA", "seqCG" };
        const sequences_type seqs2{ "AAAAAAAAAAAAAAAAAAAA"_dna4, "CGCGCGCGCGCGCGCGCGCG"_dna4 };

        // Write additional IBF
        std::string        base_prefix2{ folder_prefix + "base_build2" };
        GanonBuild::Config cfg_build;
        cfg_build.bin_size_bits      = 5000;
        cfg_build.quiet              = true;
        cfg_build.kmer_size          = 10;
        cfg_build.output_filter_file = base_prefix2 + ".ibf";
        cfg_build.reference_files    = aux::write_sequences_files( base_prefix2, "fasta", seqs2, ids2 );
        REQUIRE( GanonBuild::run( cfg_build ) );

        // Write map for base filter and additional filter
        // A2 for sequence of As in the 2nd filter
        config_classify::write_map( base_prefix + ".map", { { "0", "A" }, { "1", "C" }, { "2", "T" }, { "3", "G" } } );
        config_classify::write_map( base_prefix2 + ".map", { { "0", "A2" }, { "1", "CG" } } );

        // Write tax for base filter and  additional filter
        //        1
        //      ATCG
        //    AT    CG
        // A2 A T  C  G
        config_classify::write_tax( base_prefix + ".tax",
                                    { { "A", "AT" },
                                      { "C", "CG" },
                                      { "T", "AT" },
                                      { "G", "CG" },
                                      { "CG", "ATCG" },
                                      { "AT", "ATCG" },
                                      { "ATCG", "1" } } );
        config_classify::write_tax( base_prefix2 + ".tax",
                                    { { "A2", "AT" }, { "CG", "ATCG" }, { "AT", "ATCG" }, { "ATCG", "1" } } );

        SECTION( "with one level" )
        {
            std::string prefix{ folder_prefix + "multiple_ibf_wo_hiearchy" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix + ".ibf", base_prefix2 + ".ibf" };
            cfg.single_reads = { folder_prefix + "rA.fasta", folder_prefix + "rCG.fasta" };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Report targets as hiearchy label + filter id + bin id
            REQUIRE( res.all["readA"].size() == 3 );
            REQUIRE( res.all["readA"]["H1-0-0"] == 5 ); // A
            REQUIRE( res.all["readA"]["H1-0-2"] == 5 ); // T (rev.comp.)
            REQUIRE( res.all["readA"]["H1-1-0"] == 5 ); // A second filter
            REQUIRE( res.all["readCG"].size() == 1 );
            REQUIRE( res.all["readCG"]["H1-1-1"] == 5 ); // CG second filter

            SECTION( "--map" )
            {
                std::string prefix{ folder_prefix + "multiple_ibf_wo_hiearchy_map" };
                cfg.output_prefix = prefix;
                cfg.map           = { base_prefix + ".map", base_prefix2 + ".map" };

                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                // Matches on target (max. k-mer count)
                REQUIRE( res.all["readA"].size() == 3 );
                REQUIRE( res.all["readA"]["A"] == 5 );
                REQUIRE( res.all["readA"]["T"] == 5 );
                REQUIRE( res.all["readA"]["A2"] == 5 );
                REQUIRE( res.all["readCG"].size() == 1 );
                REQUIRE( res.all["readCG"]["CG"] == 5 );

                SECTION( "--tax" )
                {
                    std::string prefix{ folder_prefix + "multiple_ibf_wo_hiearchy_map_tax" };
                    cfg.output_prefix = prefix;
                    cfg.tax           = { base_prefix + ".tax", base_prefix2 + ".tax" };

                    REQUIRE( GanonClassify::run( cfg ) );
                    config_classify::Res res{ cfg };
                    config_classify::sanity_check( cfg, res );

                    // All matches on targets from map
                    REQUIRE( res.all["readA"].size() == 3 );
                    REQUIRE( res.all["readA"]["A"] == 5 );
                    REQUIRE( res.all["readA"]["T"] == 5 );
                    REQUIRE( res.all["readA"]["A2"] == 5 );
                    REQUIRE( res.all["readCG"].size() == 1 );
                    REQUIRE( res.all["readCG"]["CG"] == 5 );

                    // LCA matches from tax
                    REQUIRE( res.lca["readA"].size() == 1 );
                    REQUIRE( res.lca["readA"]["AT"] == 5 );
                    REQUIRE( res.lca["readCG"].size() == 1 );
                    REQUIRE( res.lca["readCG"]["CG"] == 5 );
                }
            }
        }

        SECTION( "with two levels" )
        {
            std::string prefix{ folder_prefix + "multiple_ibf_w_hiearchy" };
            auto        cfg      = config_classify::defaultConfig( prefix );
            cfg.ibf              = { base_prefix + ".ibf", base_prefix2 + ".ibf" };
            cfg.single_reads     = { folder_prefix + "rA.fasta", folder_prefix + "rCG.fasta" };
            cfg.hierarchy_labels = { "one", "two" };
            cfg.output_single    = true;

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Report targets as hiearchy label + filter id + bin id
            // same as paired and single C and G
            REQUIRE( res.all["readA"].size() == 2 );
            REQUIRE( res.all["readA"]["one-0-2"] == 5 );
            REQUIRE( res.all["readA"]["one-0-0"] == 5 );
            REQUIRE( res.all["readCG"].size() == 1 );
            REQUIRE( res.all["readCG"]["two-0-1"] == 5 );

            SECTION( "--map" )
            {
                std::string prefix{ folder_prefix + "multiple_ibf_w_hiearchy_map" };
                cfg.output_prefix = prefix;
                cfg.map           = { base_prefix + ".map", base_prefix2 + ".map" };

                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                // Matches on target (max. k-mer count)
                REQUIRE( res.all["readA"].size() == 2 );
                REQUIRE( res.all["readA"]["A"] == 5 );
                REQUIRE( res.all["readA"]["T"] == 5 );
                REQUIRE( res.all["readA"]["A2"] == 0 ); // Do not match A2 (second hiearchy, readA already classified)
                REQUIRE( res.all["readCG"].size() == 1 );
                REQUIRE( res.all["readCG"]["CG"] == 5 );

                SECTION( "--tax" )
                {
                    std::string prefix{ folder_prefix + "multiple_ibf_w_hiearchy_map_tax" };
                    cfg.output_prefix = prefix;
                    cfg.tax           = { base_prefix + ".tax", base_prefix2 + ".tax" };

                    REQUIRE( GanonClassify::run( cfg ) );
                    config_classify::Res res{ cfg };
                    config_classify::sanity_check( cfg, res );

                    // All matches on targets from map
                    REQUIRE( res.all["readA"].size() == 2 );
                    REQUIRE( res.all["readA"]["A"] == 5 );
                    REQUIRE( res.all["readA"]["T"] == 5 );
                    REQUIRE( res.all["readA"]["A2"]
                             == 0 ); // Do not match A2 (second hiearchy, readA already classified)
                    REQUIRE( res.all["readCG"].size() == 1 );
                    REQUIRE( res.all["readCG"]["CG"] == 5 );

                    // LCA matches from tax
                    REQUIRE( res.lca["readA"].size() == 1 );
                    REQUIRE( res.lca["readA"]["AT"] == 5 );
                    REQUIRE( res.lca["readCG"].size() == 1 );
                    REQUIRE( res.lca["readCG"]["CG"] == 5 );
                }
            }

            SECTION( "without --output-single" )
            {
                std::string prefix{ folder_prefix + "multiple_ibf_wo_output_single" };
                auto        cfg      = config_classify::defaultConfig( prefix );
                cfg.ibf              = { base_prefix + ".ibf", base_prefix2 + ".ibf" };
                cfg.single_reads     = { folder_prefix + "rA.fasta", folder_prefix + "rCG.fasta" };
                cfg.hierarchy_labels = { "one", "two" };
                cfg.output_single    = false;

                REQUIRE( GanonClassify::run( cfg ) );
                REQUIRE( std::filesystem::file_size( prefix + ".one.all" ) > 0 );
                REQUIRE( std::filesystem::file_size( prefix + ".two.all" ) > 0 );
            }
        }
    }
}

SCENARIO( "classifying reads with errors", "[ganon-classify][with-errors]" )
{

    std::string folder_prefix{ "ganon-classify-w-errors/" };
    std::filesystem::create_directory( folder_prefix );

    // Reads (12bp)
    aux::write_sequences( folder_prefix + "rF.fasta", { "CTCGTGTTTCCT"_dna4 }, { "readF" } );
    // RevCom from GGGCCTCTTGGT
    aux::write_sequences( folder_prefix + "rR.fasta", { "ACCAAGAGGCCC"_dna4 }, { "readR" } );

    /*
    all unknow chars (not ACTG) are replace by A, causing some unexpected assignments

    CTCGTGTTTCCT----GGGCCTCTTGGT  e0
    CTC-TGTTTCCT----GGGCCTCTTGGT  e1F
    CTC-TGTTTCCT----GGG-CTCTTGGT  e1F_e1R
    CTC-TGTTTCCT----GGG-CTCT-GGT  e1F_e2R
    CTC-TGTT-CCT----GGG-CTCTTGGT  e2F_e1R
    CTC-TGTT-CCT----GGG-CTCT-GGT  e2F_e2R

          max 4-mer count
    F  e0 9
    F  e1 5
    F  e2 1
    FR e0 18
    FR e1 14
    FR e2 10
    FR e3 6
    FR e4 2
    FR e5 2
    */
    // Sequences (28bp), error rates on references based on k = 4
    const ids_type       ids{ "e0", "e1F", "e1F_e1R", "e1F_e2R", "e2F_e1R", "e2F_e2R", "e3F_e3R" };
    const sequences_type seqs{ "CTCGTGTTTCCT----GGGCCTCTTGGT"_dna4, "CTC-TGTTTCCT----GGGCCTCTTGGT"_dna4,
                               "CTC-TGTTTCCT----GGG-CTCTTGGT"_dna4, "CTC-TGTTTCCT----GGG-CTCT-GGT"_dna4,
                               "CTC-TGTT-CCT----GGG-CTCTTGGT"_dna4, "CTC-TGTT-CCT----GGG-CTCT-GGT"_dna4 };

    std::string        base_prefix{ folder_prefix + "base_build3" };
    GanonBuild::Config cfg_build;
    cfg_build.bin_size_bits      = 5000;
    cfg_build.quiet              = true;
    cfg_build.kmer_size          = 4;
    cfg_build.output_filter_file = base_prefix + ".ibf";
    cfg_build.reference_files    = aux::write_sequences_files( base_prefix, "fasta", seqs, ids );
    REQUIRE( GanonBuild::run( cfg_build ) );

    // write map for base filter
    config_classify::write_map( base_prefix + ".map",
                                { { "0", "e0" },
                                  { "1", "e1F" },
                                  { "2", "e1F_e1R" },
                                  { "3", "e1F_e2R" },
                                  { "4", "e2F_e1R" },
                                  { "5", "e2F_e2R" } } );


    SECTION( "--abs-cutoff 0 --abs-filter 0" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_0_abs_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 0 };
        cfg.abs_filter   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Should match only e0
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0
            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
        }
    }

    SECTION( "--abs-cutoff 1 --abs-filter 0" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_1_abs_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 1 };
        cfg.abs_filter   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Should match only e0 due to strata filter
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0 due to strata filter
            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
        }
    }

    SECTION( "--abs-cutoff 2 --abs-filter 1" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_2_abs_filter_1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 2 };
        cfg.abs_filter   = { 1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Match without errors + 1 (strata)
        REQUIRE( res.all["readF"].size() == 4 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Match first and second best
            REQUIRE( res.all["readF"].size() == 2 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
        }
    }

    SECTION( "--abs-cutoff -1 (OFF) --abs-filter 2" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_-1_abs_filter_2" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { -1 };
        cfg.abs_filter   = { 2 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Match all, keep top 2 (strata)
        REQUIRE( res.all["readF"].size() == 6 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        REQUIRE( res.all["readF"]["e2F_e1R"] == 1 );
        REQUIRE( res.all["readF"]["e2F_e2R"] == 1 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Match first and second best
            REQUIRE( res.all["readF"].size() == 3 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
        }
    }

    SECTION( "--abs-cutoff 4 --abs-filter -1" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_4_abs_filter_-1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 4 };
        cfg.abs_filter   = { -1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches up-to 3 errors, without strata (-1)
        REQUIRE( res.all["readF"].size() == 6 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        REQUIRE( res.all["readF"]["e2F_e1R"] == 1 );
        REQUIRE( res.all["readF"]["e2F_e2R"] == 1 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0, not enough mers to match second best (14)
            REQUIRE( res.all["readF"].size() == 6 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e1R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e2R"] == 2 );
        }
    }

    SECTION( "--abs-cutoff -1 (OFF) --abs-filter -1 (OFF)" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_-1_abs_filter_-1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { -1 };
        cfg.abs_filter   = { -1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Match all, keep top 2 (strata) which is everything
        REQUIRE( res.all["readF"].size() == 6 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        REQUIRE( res.all["readF"]["e2F_e1R"] == 1 );
        REQUIRE( res.all["readF"]["e2F_e2R"] == 1 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 6 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e1R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e2R"] == 2 );
        }
    }

    SECTION( "--rel-cutoff 0.2 --abs-filter 0" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.2_abs_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0.2 };
        cfg.abs_filter   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Should match only e0 due to strata filter
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0 due to strata filter
            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
        }
    }

    SECTION( "--rel-cutoff 1 --abs-filter 0" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_1_abs_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 1 };
        cfg.abs_filter   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Should match only e0 matching all k-mers
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0 due to strata filter
            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
        }
    }

    SECTION( "--rel-cutoff 0.5 --abs-filter -1" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.5_abs_filter_-1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0.5 }; // 50%
        cfg.abs_filter   = { -1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Here there's a bigger difference between single and paired reads
        // due to percentage of full read being covered by kmers
        // single read matches only with one error, while paired matches up-to 2 errors
        REQUIRE( res.all["readF"].size() == 4 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0 due to strata filter
            REQUIRE( res.all["readF"].size() == 3 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
        }
    }

    SECTION( "--rel-cutoff 0 (OFF) --abs-filter -1" )
    {
        // should output every k-mer match for any k-mer count
        std::string prefix{ folder_prefix + "rel_cutoff_0_abs_filter_-1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0 };
        cfg.abs_filter   = { -1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Here there's a bigger difference between single and paired reads
        // due to percentage of full read being covered by kmers
        // single read matches only with one error, while paired matches up-to 2 errors
        REQUIRE( res.all["readF"].size() == 6 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        REQUIRE( res.all["readF"]["e2F_e1R"] == 1 );
        REQUIRE( res.all["readF"]["e2F_e2R"] == 1 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0 due to strata filter
            REQUIRE( res.all["readF"].size() == 6 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e1R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e2R"] == 2 );
        }
    }

    SECTION( "--rel-cutoff 0 (OFF) --abs-filter -1 (OFF)" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0_abs_filter_-1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0 };
        cfg.abs_filter   = { -1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 6 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        REQUIRE( res.all["readF"]["e2F_e1R"] == 1 );
        REQUIRE( res.all["readF"]["e2F_e2R"] == 1 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 6 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e1R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e2R"] == 2 );
        }
    }

    SECTION( "--abs-cutoff 2 --rel-filter 0" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_2_rel_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 2 };
        cfg.rel_filter   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
        }
    }

    SECTION( "--abs-cutoff 3 --rel-filter 0.5" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_3_rel_filter_0.5" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 3 };
        cfg.rel_filter   = { 0.5 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 4 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );


        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 3 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
        }
    }

    SECTION( "--abs-cutoff 1 --rel-filter 1 (OFF)" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_1_rel_filter_1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 1 };
        cfg.rel_filter   = { 1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 4 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 2 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
        }
    }

    SECTION( "--abs-cutoff -1 (OFF) --rel-filter 0.2" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_-1_rel_filter_0.2" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { -1 };
        cfg.rel_filter   = { 0.2 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // 9 - (ceil(9*0.2)) = 7 threshold
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // 9 - (ceil(9*0.2)) + 9 - (ceil(9*0.2)) = 14 threshold
            REQUIRE( res.all["readF"].size() == 2 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
        }
    }

    SECTION( "--abs-cutoff -1 (OFF) --rel-filter 1 (OFF)" )
    {
        std::string prefix{ folder_prefix + "abs_cutoff_-1_rel_filter_1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { -1 };
        cfg.rel_filter   = { 1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 6 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        REQUIRE( res.all["readF"]["e2F_e1R"] == 1 );
        REQUIRE( res.all["readF"]["e2F_e2R"] == 1 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 6 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e1R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e2R"] == 2 );
        }
    }

    SECTION( "--rel-cutoff 0.45 --rel-filter 0" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.45_rel_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0.45 };
        cfg.rel_filter   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // ceil(9*0.45) = 5 threhsold, 0 strata
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
        }
    }

    SECTION( "--rel-cutoff 0.2 --rel-filter 0.5" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.2_rel_filter_0.5" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0.2 };
        cfg.rel_filter   = { 0.5 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // ceil(9*0.2) = 2 threhsold cutoff
        // 9 - ceil(9*0.5) = 4 threshold filter
        REQUIRE( res.all["readF"].size() == 4 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 3 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
        }
    }

    SECTION( "--rel-cutoff 0.6 --rel-filter 1 (OFF)" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.6_rel_filter_1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0.6 };
        cfg.rel_filter   = { 1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // ceil(9*0.6) = 6 threhsold cutoff
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 2 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
        }
    }

    SECTION( "--rel-cutoff 0 (OFF) --rel-filter 0.3" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_1_rel_filter_0.3" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 0.3 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // 9 - ceil(9*0.3) = 6 threhsold cutoff
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 2 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
        }
    }

    SECTION( "--rel-cutoff 0 (OFF) --rel-filter 1 (OFF)" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_1_rel_filter_1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rF.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 1 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 6 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e1R"] == 5 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );
        REQUIRE( res.all["readF"]["e2F_e1R"] == 1 );
        REQUIRE( res.all["readF"]["e2F_e2R"] == 1 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 6 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e1R"] == 6 );
            REQUIRE( res.all["readF"]["e2F_e2R"] == 2 );
        }
    }

    SECTION( "--offset" )
    {
        /*
              offset=1 offset=2 offset=3 offset=4
        F  e0 9        5        4        3
        F  e1 5        3        3        2
        F  e2 1        1        2        1

        FR e0 18       10       8        6
        FR e1 14       8        7        5
        FR e2 10       6        6        4
        FR e3 6        4        4        3
        FR e4 2        2        3        2
        FR e5 2        2        2        1
        */
        SECTION( "--offset 2 --abs-cutoff 0 --abs-filter 0" )
        {
            std::string prefix{ folder_prefix + "offset_2_abs_cutoff_0_abs_filter_0" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.map          = { base_prefix + ".map" };
            cfg.single_reads = { folder_prefix + "rF.fasta" };
            cfg.kmer_size    = { 4 };
            cfg.abs_cutoff   = { 0 };
            cfg.abs_filter   = { 0 };
            cfg.offset       = { 2 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0
            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 5 );

            SECTION( "--paired-reads" )
            {
                prefix            = prefix + "_paired";
                cfg.output_prefix = prefix;
                cfg.single_reads  = {};
                cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };
                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                // Should match only e0
                REQUIRE( res.all["readF"].size() == 1 );
                REQUIRE( res.all["readF"]["e0"] == 10 );
            }
        }

        SECTION( "--offset 3 --abs-cutoff 4 --abs-filter -1" )
        {
            std::string prefix{ folder_prefix + "offset_3_abs_cutoff_4_abs_filter_-1" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.map          = { base_prefix + ".map" };
            cfg.single_reads = { folder_prefix + "rF.fasta" };
            cfg.kmer_size    = { 4 };
            cfg.abs_cutoff   = { 4 };
            cfg.offset       = { 3 };
            cfg.abs_filter   = { -1 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0
            REQUIRE( res.all["readF"].size() == 4 );
            REQUIRE( res.all["readF"]["e0"] == 4 );
            REQUIRE( res.all["readF"]["e1F"] == 2 );     // lost first 4-mer to offset
            REQUIRE( res.all["readF"]["e1F_e1R"] == 2 ); // lost first 4-mer to offset
            REQUIRE( res.all["readF"]["e1F_e2R"] == 2 ); // lost first 4-mer to offset

            SECTION( "--paired-reads" )
            {
                prefix            = prefix + "_paired";
                cfg.output_prefix = prefix;
                cfg.single_reads  = {};
                cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };
                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                // Should match only e0
                REQUIRE( res.all["readF"].size() == 3 );
                REQUIRE( res.all["readF"]["e0"] == 8 );
                REQUIRE( res.all["readF"]["e1F"] == 6 );     // lost first 4-mer to offset
                REQUIRE( res.all["readF"]["e1F_e1R"] == 4 ); // lost 2 4-mers to offset
            }
        }

        SECTION( "--offset 3 --rel-cutoff 0.3 --abs-filter 0" )
        {
            std::string prefix{ folder_prefix + "offset_3_rel_cutoff_0.3_abs_filter_0" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.map          = { base_prefix + ".map" };
            cfg.single_reads = { folder_prefix + "rF.fasta" };
            cfg.kmer_size    = { 4 };
            cfg.rel_cutoff   = { 0.3 };
            cfg.abs_filter   = { 0 };
            cfg.offset       = { 3 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0
            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 4 );

            SECTION( "--paired-reads" )
            {
                prefix            = prefix + "_paired";
                cfg.output_prefix = prefix;
                cfg.single_reads  = {};
                cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };
                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                // Should match only e0
                REQUIRE( res.all["readF"].size() == 1 );
                REQUIRE( res.all["readF"]["e0"] == 8 );
            }
        }

        SECTION( "--offset 4 --rel-cutoff 0.5 --rel-filter 0.5" )
        {
            std::string prefix{ folder_prefix + "offset_4_rel_cutoff_0.5_rel_filter_0.5" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.map          = { base_prefix + ".map" };
            cfg.single_reads = { folder_prefix + "rF.fasta" };
            cfg.kmer_size    = { 4 };
            cfg.rel_cutoff   = { 0.5 };
            cfg.rel_filter   = { 0.5 };
            cfg.offset       = { 4 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 4 );
            REQUIRE( res.all["readF"]["e0"] == 3 );
            REQUIRE( res.all["readF"]["e1F"] == 2 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 2 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 2 );

            SECTION( "--paired-reads" )
            {
                prefix            = prefix + "_paired";
                cfg.output_prefix = prefix;
                cfg.single_reads  = {};
                cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };
                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                REQUIRE( res.all["readF"].size() == 5 );
                REQUIRE( res.all["readF"]["e0"] == 6 );
                REQUIRE( res.all["readF"]["e1F"] == 5 );
                REQUIRE( res.all["readF"]["e1F_e1R"] == 4 );
                REQUIRE( res.all["readF"]["e1F_e2R"] == 3 );
                REQUIRE( res.all["readF"]["e2F_e1R"] == 3 );
            }
        }
    }


    SECTION( "--window-size" )
    {
        // build with --window-size
        std::string        base_prefix_ws{ folder_prefix + "base_build_window_size" };
        GanonBuild::Config cfg_build_ws;
        cfg_build_ws.bin_size_bits      = 5000;
        cfg_build_ws.quiet              = true;
        cfg_build_ws.kmer_size          = 4;
        cfg_build_ws.window_size        = 6;
        cfg_build_ws.output_filter_file = base_prefix_ws + ".ibf";
        cfg_build_ws.reference_files    = aux::write_sequences_files( base_prefix_ws, "fasta", seqs, ids );
        REQUIRE( GanonBuild::run( cfg_build_ws ) );

        SECTION( "--window-size 6 --rel-cutoff 1 --rel-filter 0" )
        {
            std::string prefix{ folder_prefix + "window_size_6_rel_cutoff_1_rel_filter_0" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix_ws + ".ibf" };
            cfg.map          = { base_prefix + ".map" };
            cfg.single_reads = { folder_prefix + "rF.fasta" };
            cfg.kmer_size    = { 4 };
            cfg.rel_cutoff   = { 1 };
            cfg.rel_filter   = { 0 };
            cfg.window_size  = { 6 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0
            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 4 );

            SECTION( "--paired-reads" )
            {
                prefix            = prefix + "_paired";
                cfg.output_prefix = prefix;
                cfg.single_reads  = {};
                cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };
                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                // Should match only e0
                REQUIRE( res.all["readF"].size() == 1 );
                REQUIRE( res.all["readF"]["e0"] == 8 );
            }
        }

        SECTION( "--window-size 6 --rel-cutoff 0 --rel-filter 1" )
        {
            std::string prefix{ folder_prefix + "window_size_6_rel_cutoff_0_rel_filter_1" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix_ws + ".ibf" };
            cfg.map          = { base_prefix + ".map" };
            cfg.single_reads = { folder_prefix + "rF.fasta" };
            cfg.kmer_size    = { 4 };
            cfg.rel_cutoff   = { 0 };
            cfg.rel_filter   = { 1 };
            cfg.window_size  = { 6 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 4 );
            REQUIRE( res.all["readF"]["e0"] == 4 );
            REQUIRE( res.all["readF"]["e1F"] == 2 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 2 );
            REQUIRE( res.all["readF"]["e1F_e2R"] == 2 );

            SECTION( "--paired-reads" )
            {
                prefix            = prefix + "_paired";
                cfg.output_prefix = prefix;
                cfg.single_reads  = {};
                cfg.paired_reads  = { folder_prefix + "rF.fasta", folder_prefix + "rR.fasta" };
                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                REQUIRE( res.all["readF"].size() == 5 );
                REQUIRE( res.all["readF"]["e0"] == 8 );
                REQUIRE( res.all["readF"]["e1F"] == 6 );
                REQUIRE( res.all["readF"]["e1F_e1R"] == 4 );
                REQUIRE( res.all["readF"]["e1F_e2R"] == 2 );
                REQUIRE( res.all["readF"]["e2F_e1R"] == 2 );
            }
        }
    }

    SECTION( "reads with errors and --abs-cutoff 1 --abs-filter 0" )
    {


        // using same filter twice
        // reads with 1 error at beginning
        aux::write_sequences( folder_prefix + "rFe1.fasta", { "CTCGTGTTTCC-"_dna4 }, { "readFe1" } );
        aux::write_sequences( folder_prefix + "rRe1.fasta", { "ACCAAGAGGCC-"_dna4 }, { "readRe1" } );

        std::string prefix{ folder_prefix + "reads_with_error_abs_cutoff_1" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.map          = { base_prefix + ".map" };
        cfg.single_reads = { folder_prefix + "rFe1.fasta" };
        cfg.kmer_size    = { 4 };
        cfg.abs_cutoff   = { 1 };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readFe1"].size() == 1 );
        REQUIRE( res.all["readFe1"]["e0"] == 8 );


        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "rFe1.fasta", folder_prefix + "rRe1.fasta" };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0
            REQUIRE( res.all["readFe1"].size() == 1 );
            REQUIRE( res.all["readFe1"]["e0"] == 16 );
        }
    }
}