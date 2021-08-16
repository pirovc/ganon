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
typedef std::map< std::string,  std::string >                      TTax;

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


SCENARIO( "classifying reads without errors", "[ganon-classify]" )
{

    // Reads (14bp)
    aux::write_sequences( "rA.fasta", { "AAAAAAAAAAAAAA"_dna5 }, { "readA" } );
    aux::write_sequences( "rC.fasta", { "CCCCCCCCCCCCCC"_dna5 }, { "readC" } );
    aux::write_sequences( "rT.fasta", { "TTTTTTTTTTTTTT"_dna5 }, { "readT" } );
    aux::write_sequences( "rG.fasta", { "GGGGGGGGGGGGGG"_dna5 }, { "readG" } );

    // Sequences (20bp)
    const ids_type       ids{ "seqA", "seqC", "seqT", "seqG" };
    const sequences_type seqs{ "AAAAAAAAAAAAAAAAAAAA"_dna5,
                               "CCCCCCCCCCCCCCCCCCCC"_dna5,
                               "TTTTTTTTTTTTTTTTTTTT"_dna5,
                               "GGGGGGGGGGGGGGGGGGGG"_dna5 };


    std::string        base_prefix{ "classify_base_build" };
    GanonBuild::Config cfg_build;
    cfg_build.bin_size_bits      = 5000;
    cfg_build.quiet              = true;
    cfg_build.kmer_size          = 10;
    cfg_build.output_filter_file = base_prefix + ".ibf";
    cfg_build.reference_files    = aux::write_sequences_files( base_prefix, "fasta", seqs, ids );
    REQUIRE( GanonBuild::run( cfg_build ) );

    SECTION( "with default config." )
    {
        std::string prefix{ "classify_base_build" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { "rA.fasta" };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches on binids (Rev. Comp.)
        REQUIRE( res.all["readA"]["0"] == 5 );
        REQUIRE( res.all["readA"]["1"] == 0 );
        REQUIRE( res.all["readA"]["2"] == 5 );
        REQUIRE( res.all["readA"]["3"] == 0 );
    }

    SECTION( "with --map" )
    {
        std::string prefix{ "classify_map" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { "rA.fasta" };
        config_classify::write_map( "classify_map.map",
                                    { { "0", "AorT" }, { "1", "CorG" }, { "2", "AorT" }, { "3", "CorG" } } );
        cfg.map = { "classify_map.map" };

        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // Matches on target (max. k-mer count)
        REQUIRE( res.all["readA"]["AorT"] == 5 );
        REQUIRE( res.all["readA"]["CorG"] == 0 );
    }

    SECTION( "with --map and --tax" )
    {
        std::string prefix{ "classify_map_tax" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { "rA.fasta" };
        config_classify::write_map( "classify_map_tax.map",
                                    { { "0", "A" }, { "1", "C" }, { "2", "T" }, { "3", "G" } } );
        cfg.map = { "classify_map_tax.map" };

        //     1
        //    ATCG
        //  AT    CG
        // A  T  C  G
        config_classify::write_tax( "classify_map_tax.tax",
                                    { { "A", "AT" }, { "C", "CG" }, { "T", "AT" }, { "G", "CG" },
                                      { "CG", "ATCG" },  { "AT", "ATCG" }, { "ATCG", "1" } } );
        cfg.tax = { "classify_map_tax.tax" };
        
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // All matches on targets from map
        REQUIRE( res.all["readA"]["A"] == 5 );
        REQUIRE( res.all["readA"]["T"] == 5 );

        // LCA matches from tax
        REQUIRE( res.lca["readA"]["AT"] == 5 );
    }

    SECTION( "with --paired-reads" )
    {
        std::string prefix{ "classify_paired_reads" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.paired_reads = { "rA.fasta", "rT.fasta" };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // will report header of first pair "readA" and match the rev.comp. of readT (and opposite)
        REQUIRE( res.all["readA"]["0"] == 10 );
        REQUIRE( res.all["readA"]["1"] == 0 );
        REQUIRE( res.all["readA"]["2"] == 10 );
        REQUIRE( res.all["readA"]["3"] == 0 );
    }

    SECTION( "with --single-reads and --paired-reads" )
    {
        std::string prefix{ "classify_single_paired_reads" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { "rC.fasta", "rG.fasta" };
        cfg.paired_reads = { "rA.fasta", "rT.fasta" };
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{ cfg };
        config_classify::sanity_check( cfg, res );

        // same as paired and single C and G
        REQUIRE( res.all["readA"]["0"] == 10 );
        REQUIRE( res.all["readC"]["1"] == 5 );
        REQUIRE( res.all["readG"]["1"] == 5 );
        REQUIRE( res.all["readA"]["2"] == 10 );
        REQUIRE( res.all["readC"]["3"] == 5 );
        REQUIRE( res.all["readG"]["3"] == 5 );
    }

    SECTION( "with wrong --kmer-size" )
    {
        std::string prefix{ "classify_wrong_kmer_size" };
        auto        cfg  = config_classify::defaultConfig( prefix );
        cfg.ibf          = { base_prefix + ".ibf" };
        cfg.single_reads = { "rA.fasta" };
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
        std::string prefix{ "classify_wo_output_prefix" };
        auto        cfg   = config_classify::defaultConfig( prefix );
        cfg.ibf           = { base_prefix + ".ibf" };
        cfg.single_reads  = { "rA.fasta" };
        cfg.output_prefix = "";
        REQUIRE( GanonClassify::run( cfg ) );
        // No files created
        REQUIRE_FALSE( std::filesystem::exists( "classify_wo_output_prefix.rep" ) );
        REQUIRE_FALSE( std::filesystem::exists( "classify_wo_output_prefix.lca" ) );
        REQUIRE_FALSE( std::filesystem::exists( "classify_wo_output_prefix.all" ) );
        REQUIRE_FALSE( std::filesystem::exists( "classify_wo_output_prefix.unc" ) );
    }


    SECTION( "with multiple --ibf" )
    {

        // Sequences (20bp)
        const ids_type       ids2{ "seqA", "seqN" };
        const sequences_type seqs2{ "AAAAAAAAAAAAAAAAAAAA"_dna5, "NNNNNNNNNNNNNNNNNNNN"_dna5 };

        std::string        base_prefix2{ "classify_base_build2" };
        GanonBuild::Config cfg_build;
        cfg_build.bin_size_bits      = 5000;
        cfg_build.quiet              = true;
        cfg_build.kmer_size          = 10;
        cfg_build.output_filter_file = base_prefix2 + ".ibf";
        cfg_build.reference_files    = aux::write_sequences_files( base_prefix2, "fasta", seqs, ids );
        REQUIRE( GanonBuild::run( cfg_build ) );

        SECTION( "without hiearchy" )
        {
            std::string prefix{ "classify_multiple_ibf_wo_hiearchy" };
            auto        cfg  = config_classify::defaultConfig( prefix );
            cfg.ibf          = { base_prefix + ".ibf", base_prefix2 + ".ibf" };
            cfg.single_reads = { "rA.fasta" };

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Report targets as hiearchy label + filter id + bin id
            // same as paired and single C and G
            REQUIRE( res.all["readA"]["1_default-1-2"] == 5 );
            REQUIRE( res.all["readA"]["1_default-1-0"] == 5 );
            REQUIRE( res.all["readA"]["1_default-0-0"] == 5 );
            REQUIRE( res.all["readA"]["1_default-0-2"] == 5 );
        }

        SECTION( "with hiearchy" )
        {
            std::string prefix{ "classify_multiple_ibf_w_hiearchy" };
            auto        cfg      = config_classify::defaultConfig( prefix );
            cfg.ibf              = { base_prefix + ".ibf", base_prefix2 + ".ibf" };
            cfg.single_reads     = { "rA.fasta" };
            cfg.hierarchy_labels = { "one", "two" };
            cfg.output_single    = true;

            REQUIRE( GanonClassify::run( cfg ) );
            config_classify::Res res{ cfg };
            config_classify::sanity_check( cfg, res );

            // Report targets as hiearchy label + filter id + bin id
            // same as paired and single C and G
            REQUIRE( res.all["readA"]["one-0-2"] == 5 );
            REQUIRE( res.all["readA"]["one-0-0"] == 5 );
        }
    }
}

SCENARIO( "classifying reads with errors", "[ganon-classify]" )
{
    SECTION( "..." )
    {
        REQUIRE( true );
    }
}