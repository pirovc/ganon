#include "aux/Aux.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
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
    cfg.verbose             = false;
    cfg.quiet               = true;
    return cfg;
}


typedef std::map< std::string, std::map< std::string, uint32_t > > TOut;
typedef std::vector< std::string >                                 TUnc;
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
            //   0 hierarchy_label <tab>
            //   1 target <tab>
            //   2 matches <tab>
            //   3 unique reads <tab>
            //   4 lca reads <tab>
            //   5 [rank <tab>]
            //   6 [name]
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
            lca = parse_all_lca( cfg.output_prefix + ".one", lines_lca );
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
    aux::sequences_type reads{
        "AAAAAAAAAAAAAA"_dna4, "CCCCCCCCCCCCCC"_dna4, "TTTTTTTTTTTTTT"_dna4, "GGGGGGGGGGGGGG"_dna4
    };

    auto seqtarget_reads = aux::SeqTarget( folder_prefix, reads, {}, { "readA", "readC", "readT", "readG" } );
    seqtarget_reads.write_sequences_files();


    // Reference sequences (20bp)
    aux::sequences_type refs{ "AAAAAAAAAAAAAAAAAAAA"_dna4,
                              "CCCCCCCCCCCCCCCCCCCC"_dna4,
                              "TTTTTTTTTTTTTTTTTTTT"_dna4,
                              "GGGGGGGGGGGGGGGGGGGG"_dna4 };

    auto seqtarget_refs =
        aux::SeqTarget( folder_prefix, refs, { "A", "C", "T", "G" }, { "seqA", "seqC", "seqT", "seqG" } );

    std::string base_prefix{ folder_prefix + "base_build1" };
    seqtarget_refs.write_input_file( base_prefix + ".tsv" );
    seqtarget_refs.write_sequences_files();

    GanonBuild::Config cfg_build;
    cfg_build.input_file  = base_prefix + ".tsv";
    cfg_build.output_file = base_prefix + ".ibf";
    cfg_build.max_fp      = 0.01;
    cfg_build.verbose     = false;
    cfg_build.quiet       = true;
    cfg_build.kmer_size   = 10;
    cfg_build.window_size = 10;
    REQUIRE( GanonBuild::run( cfg_build ) );

    SECTION( "--verbose" )
    {
        std::string prefix{ folder_prefix + "verbose" };

        // Redirect cerr to file
        std::streambuf* backup_cerr;
        std::ofstream   filestr;
        filestr.open( prefix + ".log" );
        backup_cerr = std::cerr.rdbuf();
        std::cerr.rdbuf( filestr.rdbuf() );

        auto cfg         = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readA.fasta" };
        cfg.verbose      = true;

        REQUIRE( GanonClassify::run( cfg ) );

        // restore cerr and close file
        std::cerr.rdbuf( backup_cerr );
        filestr.close();

        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // check if there's output on verbose log
        REQUIRE_FALSE( aux::fileIsEmpty( prefix + ".log" ) );
    }

    SECTION( "--single-reads" )
    {
        std::string prefix{ folder_prefix + "single_reads" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readA.fasta" };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 1 };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches on binids (Rev. Comp.)
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["A"] == 5 );
        REQUIRE( res.all["readA"]["T"] == 5 );

        SECTION( "without --output-lca" )
        {
            std::string prefix{ folder_prefix + "single_reads_wo_output_lca" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.output_lca   = false;
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.single_reads = { folder_prefix + "readA.fasta" };
            cfg.rel_cutoff   = { 0 };
            cfg.rel_filter   = { 1 };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );
            REQUIRE_FALSE( std::filesystem::exists( prefix + ".one" ) );
        }
        SECTION( "without --output-all" )
        {
            std::string prefix{ folder_prefix + "default_wo_output_all" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.output_all   = false;
            cfg.ibf          = { base_prefix + ".ibf" };
            cfg.single_reads = { folder_prefix + "readA.fasta" };
            cfg.rel_cutoff   = { 0 };
            cfg.rel_filter   = { 1 };

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
        cfg.paired_reads = { folder_prefix + "readA.fasta", folder_prefix + "readT.fasta" };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 1 };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // will report header of first pair "readA" and match the rev.comp. of readT (and opposite)
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["A"] == 10 );
        REQUIRE( res.all["readA"]["T"] == 10 );
    }

    SECTION( "--single-reads and --paired-reads" )
    {
        std::string prefix{ folder_prefix + "single_paired_reads" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readC.fasta", folder_prefix + "readG.fasta" };
        cfg.paired_reads = { folder_prefix + "readA.fasta", folder_prefix + "readT.fasta" };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 1 };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // same as paired and single C and G
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readC"].size() == 2 );
        REQUIRE( res.all["readG"].size() == 2 );
        REQUIRE( res.all["readA"]["A"] == 10 );
        REQUIRE( res.all["readC"]["C"] == 5 );
        REQUIRE( res.all["readG"]["C"] == 5 );
        REQUIRE( res.all["readA"]["T"] == 10 );
        REQUIRE( res.all["readC"]["C"] == 5 );
        REQUIRE( res.all["readG"]["C"] == 5 );
    }


    SECTION( "--tax" )
    {
        std::string prefix{ folder_prefix + "tax" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readA.fasta" };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 1 };

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


    SECTION( "incomplete --tax" )
    {
        std::string prefix{ folder_prefix + "incomplete_tax" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readA.fasta" };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 1 };

        //     1
        //    ATCG
        //  AT    CG
        // -  T  C  G
        config_classify::write_tax( prefix + ".tax",
                                    { { "C", "CG" },
                                      { "T", "AT" },
                                      { "G", "CG" },
                                      { "CG", "ATCG" },
                                      { "AT", "ATCG" },
                                      { "ATCG", "1" } } ); // missing  { "A", "AT" }
        cfg.tax = { prefix + ".tax" };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches are correctly assigned on target
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["A"] == 5 );
        REQUIRE( res.all["readA"]["T"] == 5 );

        // A is missing from tax, so it's linked to to root (1), LCA is root (1)
        REQUIRE( res.lca["readA"].size() == 1 );
        REQUIRE( res.lca["readA"]["1"] == 5 );
    }

    SECTION( "--window-size != --kmer-size" )
    {
        // build with --window-size
        std::string        base_prefix_ws{ folder_prefix + "base_build_window_size" };
        GanonBuild::Config cfg_build_ws;
        cfg_build_ws.output_file = base_prefix_ws + ".ibf";
        cfg_build_ws.input_file  = base_prefix + ".tsv";
        cfg_build_ws.max_fp      = 0.01;
        cfg_build_ws.quiet       = true;
        cfg_build_ws.kmer_size   = 10;
        cfg_build_ws.window_size = 12;
        REQUIRE( GanonBuild::run( cfg_build_ws ) );

        std::string prefix{ folder_prefix + "window_size" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix_ws + ".ibf" };
        cfg.single_reads = { folder_prefix + "readA.fasta" };
        cfg.rel_filter   = { 1 };
        cfg.rel_cutoff   = { 0 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // one kmer counted since they share same hash
        REQUIRE( res.all["readA"].size() == 2 );
        REQUIRE( res.all["readA"]["A"] == 1 );
        REQUIRE( res.all["readA"]["T"] == 1 );

        SECTION( "--paired-reads" )
        {
            std::string prefix{ folder_prefix + "window_size_paired_reads" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix_ws + ".ibf" };
            cfg.paired_reads = { folder_prefix + "readA.fasta", folder_prefix + "readT.fasta" };
            cfg.rel_filter   = { 1 };
            cfg.rel_cutoff   = { 0 };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // will report header of first pair "readA" and match the rev.comp. of readT (and opposite)
            REQUIRE( res.all["readA"].size() == 2 );
            REQUIRE( res.all["readA"]["A"] == 2 );
            REQUIRE( res.all["readA"]["T"] == 2 );
        }
    }

    SECTION( "without --output-prefix" )
    {
        // avoid printing on test results (ignoring STDOUT)
        std::streambuf* old = std::cout.rdbuf( 0 );

        std::string prefix{ folder_prefix + "wo_output_prefix" };
        auto        cfg   = config_classify::defaultConfig( prefix );
        cfg.ibf           = { base_prefix + ".ibf" };
        cfg.single_reads  = { folder_prefix + "readA.fasta" };
        cfg.output_prefix = "";

        REQUIRE( GanonClassify::run( cfg ) );
        // No files created
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.rep" ) );
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.one" ) );
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.all" ) );
        REQUIRE_FALSE( std::filesystem::exists( "wo_output_prefix.unc" ) );

        std::cout.rdbuf( old );
    }

    SECTION( "--hierarchy-labels" )
    {
        // Additional filter with repeated sequences (As) and new sequence (CG)

        // Reads (14bp)
        aux::sequences_type reads2{ "CGCGCGCGCGCGCG"_dna4 };
        auto                seqtarget_reads2 = aux::SeqTarget( folder_prefix, reads2, { "CG" }, { "readCG" } );
        seqtarget_reads2.write_sequences_files();

        // Reference sequences (20bp)
        aux::sequences_type refs2{
            "AAAAAAAAAAAAAAAAAAAA"_dna4,
            "CGCGCGCGCGCGCGCGCGCG"_dna4,
        };

        auto        seqtarget_refs2 = aux::SeqTarget( folder_prefix, refs2, { "A2", "CG" }, { "seqA2", "seqCG" } );
        std::string base_prefix2{ folder_prefix + "base_build2" };
        seqtarget_refs2.write_input_file( base_prefix2 + ".tsv" );
        seqtarget_refs2.write_sequences_files();

        // Write additional IBF
        GanonBuild::Config cfg_build;
        cfg_build.output_file = base_prefix2 + ".ibf";
        cfg_build.input_file  = base_prefix2 + ".tsv";
        cfg_build.max_fp      = 0.01;
        cfg_build.quiet       = true;
        cfg_build.kmer_size   = 10;
        cfg_build.window_size = 10;
        REQUIRE( GanonBuild::run( cfg_build ) );


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
            cfg.single_reads = { folder_prefix + "readA.fasta", folder_prefix + "readCG.fasta" };
            cfg.rel_cutoff   = { 0 };
            cfg.rel_filter   = { 1 };
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
                std::string prefix{ folder_prefix + "multiple_ibf_wo_hiearchy_tax" };
                cfg.output_prefix = prefix;
                cfg.tax           = { base_prefix + ".tax", base_prefix2 + ".tax" };
                cfg.rel_cutoff    = { 0 };
                cfg.rel_filter    = { 1 };

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

        SECTION( "with two levels" )
        {
            std::string prefix{ folder_prefix + "multiple_ibf_w_hiearchy" };
            auto        cfg      = config_classify::defaultConfig( prefix );
            cfg.ibf              = { base_prefix + ".ibf", base_prefix2 + ".ibf" };
            cfg.single_reads     = { folder_prefix + "readA.fasta", folder_prefix + "readCG.fasta" };
            cfg.hierarchy_labels = { "one", "two" };
            cfg.output_single    = true;
            cfg.rel_cutoff       = { 0 };
            cfg.rel_filter       = { 1 };

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
                cfg.rel_cutoff    = { 0 };
                cfg.rel_filter    = { 1 };

                REQUIRE( GanonClassify::run( cfg ) );
                config_classify::Res res{ cfg };
                config_classify::sanity_check( cfg, res );

                // All matches on targets from map
                REQUIRE( res.all["readA"].size() == 2 );
                REQUIRE( res.all["readA"]["A"] == 5 );
                REQUIRE( res.all["readA"]["T"] == 5 );
                REQUIRE( res.all["readA"]["A2"] == 0 ); // Do not match A2 (second hiearchy, readA already classified)
                REQUIRE( res.all["readCG"].size() == 1 );
                REQUIRE( res.all["readCG"]["CG"] == 5 );

                // LCA matches from tax
                REQUIRE( res.lca["readA"].size() == 1 );
                REQUIRE( res.lca["readA"]["AT"] == 5 );
                REQUIRE( res.lca["readCG"].size() == 1 );
                REQUIRE( res.lca["readCG"]["CG"] == 5 );
            }


            SECTION( "without --output-single" )
            {
                std::string prefix{ folder_prefix + "multiple_ibf_wo_output_single" };
                auto        cfg = config_classify::defaultConfig( prefix );
                cfg.rel_cutoff  = { 0 };
                cfg.rel_filter  = { 1 };

                cfg.ibf              = { base_prefix + ".ibf", base_prefix2 + ".ibf" };
                cfg.single_reads     = { folder_prefix + "readA.fasta", folder_prefix + "readCG.fasta" };
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

    // Reads (14bp)
    aux::sequences_type reads{
        "CTCGTGTTTCCT"_dna4,
        "ACCAAGAGGCCC"_dna4,
    };
    auto seqtarget_reads = aux::SeqTarget( folder_prefix, reads, {}, { "readF", "readR" } );
    seqtarget_reads.write_sequences_files();

    // all unknow chars (not ACTG) are replaced by A, causing some unexpected assignments
    //
    // CTCGTGTTTCCT----GGGCCTCTTGGT  e0
    // CTC-TGTTTCCT----GGGCCTCTTGGT  e1F
    // CTC-TGTTTCCT----GGG-CTCTTGGT  e1F_e1R
    // CTC-TGTTTCCT----GGG-CTCT-GGT  e1F_e2R
    // CTC-TGTT-CCT----GGG-CTCTTGGT  e2F_e1R
    // CTC-TGTT-CCT----GGG-CTCT-GGT  e2F_e2R
    //
    //      max 4-mer count
    // F  e0 9
    // F  e1 5
    // F  e2 1
    // FR e0 18
    // FR e1 14
    // FR e2 10
    // FR e3 6
    // FR e4 2
    // FR e5 2


    // Reference sequences (20bp)
    aux::sequences_type refs{ "CTCGTGTTTCCT----GGGCCTCTTGGT"_dna4, "CTC-TGTTTCCT----GGGCCTCTTGGT"_dna4,
                              "CTC-TGTTTCCT----GGG-CTCTTGGT"_dna4, "CTC-TGTTTCCT----GGG-CTCT-GGT"_dna4,
                              "CTC-TGTT-CCT----GGG-CTCTTGGT"_dna4, "CTC-TGTT-CCT----GGG-CTCT-GGT"_dna4 };
    auto                seqtarget_refs = aux::SeqTarget( folder_prefix,
                                          refs,
                                          { "e0", "e1F", "e1F_e1R", "e1F_e2R", "e2F_e1R", "e2F_e2R", "e3F_e3R" },
                                          { "e0", "e1F", "e1F_e1R", "e1F_e2R", "e2F_e1R", "e2F_e2R", "e3F_e3R" } );

    std::string base_prefix{ folder_prefix + "base_build3" };
    seqtarget_refs.write_input_file( base_prefix + ".tsv" );
    seqtarget_refs.write_sequences_files();

    GanonBuild::Config cfg_build;
    cfg_build.output_file = base_prefix + ".ibf";
    cfg_build.input_file  = base_prefix + ".tsv";
    cfg_build.max_fp      = 0.01;
    cfg_build.quiet       = true;
    cfg_build.kmer_size   = 4;
    cfg_build.window_size = 4;
    REQUIRE( GanonBuild::run( cfg_build ) );

    SECTION( "--rel-cutoff 0.45 --rel-filter 0" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.45_rel_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readF.fasta" };
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
            cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 1 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
        }
    }

    SECTION( "--rel-cutoff 0.2 --rel-filter 0.8" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.2_rel_filter_0.8" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readF.fasta" };
        cfg.rel_cutoff   = { 0.2 };
        cfg.rel_filter   = { 0.8 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // ceil(9*0.2) = 2 threhsold cutoff
        // 9 - ceil((9-5)*0.8) = 5 threshold filter
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
            cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 3 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
        }
    }

    SECTION( "--rel-cutoff 0.2 --rel-filter 0.8 --fpr-query 1e-10" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_0.2_rel_filter_0.8_fpr_query_1e-10" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readF.fasta" };
        cfg.rel_cutoff   = { 0.2 };
        cfg.rel_filter   = { 0.8 };
        cfg.fpr_query    = { 1e-10 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 2 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };

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
        cfg.single_reads = { folder_prefix + "readF.fasta" };
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
            cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };

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
        cfg.single_reads = { folder_prefix + "readF.fasta" };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 0.3 };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // 9 - ceil((9-9)*0.3) = 9 threhsold cutoff
        REQUIRE( res.all["readF"].size() == 1 );
        REQUIRE( res.all["readF"]["e0"] == 9 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };

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
        cfg.single_reads = { folder_prefix + "readF.fasta" };
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
            cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };

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

    SECTION( "--rel-cutoff 0 (OFF) --rel-filter 1 (OFF) --fpr-query 1e-10" )
    {
        std::string prefix{ folder_prefix + "rel_cutoff_1_rel_filter_1_fpr_query_1e-10" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readF.fasta" };
        cfg.rel_cutoff   = { 0 };
        cfg.rel_filter   = { 1 };
        cfg.fpr_query    = { 1e-10 };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        REQUIRE( res.all["readF"].size() == 2 );
        REQUIRE( res.all["readF"]["e0"] == 9 );
        REQUIRE( res.all["readF"]["e1F_e2R"] == 5 );

        SECTION( "--paired-reads" )
        {
            prefix            = prefix + "_paired";
            cfg.output_prefix = prefix;
            cfg.single_reads  = {};
            cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            REQUIRE( res.all["readF"].size() == 3 );
            REQUIRE( res.all["readF"]["e0"] == 18 );
            REQUIRE( res.all["readF"]["e1F"] == 14 );
            REQUIRE( res.all["readF"]["e1F_e1R"] == 10 );
        }
    }

    SECTION( "--window-size" )
    {
        // build with --window-size
        std::string        base_prefix_ws{ folder_prefix + "base_build_window_size" };
        GanonBuild::Config cfg_build_ws;
        cfg_build_ws.output_file = base_prefix_ws + ".ibf";
        cfg_build_ws.input_file  = base_prefix + ".tsv";
        cfg_build_ws.max_fp      = 0.01;
        cfg_build_ws.quiet       = true;
        cfg_build_ws.kmer_size   = 4;
        cfg_build_ws.window_size = 6;
        REQUIRE( GanonBuild::run( cfg_build_ws ) );

        SECTION( "--window-size 6 --rel-cutoff 1 --rel-filter 0" )
        {
            std::string prefix{ folder_prefix + "window_size_6_rel_cutoff_1_rel_filter_0" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix_ws + ".ibf" };
            cfg.single_reads = { folder_prefix + "readF.fasta" };
            cfg.rel_cutoff   = { 1 };
            cfg.rel_filter   = { 0 };

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
                cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };
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
            cfg.single_reads = { folder_prefix + "readF.fasta" };
            cfg.rel_cutoff   = { 0 };
            cfg.rel_filter   = { 1 };

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
                cfg.paired_reads  = { folder_prefix + "readF.fasta", folder_prefix + "readR.fasta" };
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

    SECTION( "reads with errors and --rel-cutoff 0.7 --rel-filter 0" )
    {

        // using same filter twice
        // reads with 1 error at beginning
        aux::sequences_type reads{ "CTCGTGTTTCC-"_dna4, "ACCAAGAGGCC-"_dna4 };

        auto seqtarget_reads = aux::SeqTarget( folder_prefix, reads, {}, { "readFe1", "readRe1" } );
        seqtarget_reads.write_sequences_files();

        std::string prefix{ folder_prefix + "reads_with_error_rel_cutoff_0.7_rel_filter_0" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { folder_prefix + "readFe1.fasta" };
        cfg.rel_cutoff   = { 0.7 };
        cfg.rel_filter   = { 0 };

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
            cfg.paired_reads  = { folder_prefix + "readFe1.fasta", folder_prefix + "readRe1.fasta" };
            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Should match only e0
            REQUIRE( res.all["readFe1"].size() == 1 );
            REQUIRE( res.all["readFe1"]["e0"] == 16 );
        }
    }
}
