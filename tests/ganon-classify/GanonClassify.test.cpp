#include "aux/Aux.hpp"

#include <seqan3/core/debug_stream.hpp>

#include <ganon-classify/Config.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <ganon-build/Config.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <catch2/catch.hpp>

using namespace seqan3::literals;

namespace config_classify
{

GanonClassify::Config defaultConfig( const std::string prefix  )
{
    GanonClassify::Config cfg;
    cfg.output_prefix = prefix;
    cfg.output_all          = true;
    cfg.output_lca          = true;
    cfg.output_unclassified = true;
    cfg.threads             = 4;
    cfg.kmer_size           = { 10 };
    cfg.verbose             = false;
    cfg.quiet               = true;
    return cfg;
}

typedef std::map< std::string, std::map<std::string, uint32_t> > TOut;
typedef std::vector< std::string > TUnc;

struct Res{
    
    void parse_rep(std::string file){
        std::string   line;
        std::ifstream infile{file};
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
            if(fields[0]=="#total_classified"){
                total_classified = std::stoul( fields[1] );
            }else if(fields[0]=="#total_unclassified"){
                total_unclassified = std::stoul( fields[1] );
            }else{
                matches += std::stoul( fields[2] );
                unique_reads += std::stoul( fields[3] );
                lca_reads += std::stoul( fields[4] );
            }
        }
        total_reads = total_classified + total_unclassified;
        infile.close();
    }

    TOut parse_all_lca(std::string file, uint64_t & nlines){
        TOut out;
        std::string   line;
        std::ifstream infile{file};
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

    TUnc parse_unc(std::string file){
        TUnc out;
        std::string line;
        std::ifstream infile{file};
        // Read the next line from File untill it reaches the end.
        while (std::getline(infile, line))
        {
            out.push_back(line);
        }
        return out;
    }

    Res()
    {
    }

    Res(GanonClassify::Config& cfg){
        parse_rep(cfg.output_prefix + ".rep");
        if(cfg.output_all){
            all = parse_all_lca( cfg.output_prefix + ".all", lines_all );
        }
        if(cfg.output_lca){
            lca = parse_all_lca( cfg.output_prefix + ".lca", lines_lca );
        }
        if(cfg.output_unclassified){
           unc = parse_unc( cfg.output_prefix + ".unc" );
        }
    }

    TOut all;
    TOut lca;
    TUnc unc;
    uint64_t total_classified=0;
    uint64_t total_unclassified=0;
    uint64_t total_reads=0;
    uint64_t matches=0;
    uint64_t unique_reads=0;
    uint64_t lca_reads=0;
    uint64_t lines_all = 0;
    uint64_t lines_lca = 0;
};

void sanity_check(GanonClassify::Config const& cfg, Res const& res){
    if(cfg.output_all){
        // Output as many reads and matches as reported (.rep)
        REQUIRE( res.all.size()==res.total_classified );
        REQUIRE( res.lines_all==res.matches );
    }

    if(cfg.output_lca && !cfg.tax.empty()){
        // Output as many lca reads as reported (.rep)
        REQUIRE( res.lca.size()==res.total_classified );
        REQUIRE( res.lines_lca==res.total_classified );
    }

    if(cfg.output_unclassified){
        // Output as many unclassified reads as reported (.rep)
        REQUIRE( res.unc.size()==res.total_unclassified );
    }
}

} // namespace config_classify


SCENARIO( "classifying reads without errors", "[ganon-classify]" )
{

    const ids_type       ids{ "seqA", "seqC", "seqT", "seqG" };
    const sequences_type seqs{ "AAAAAAAAAAAAAAAAAAAA"_dna5,
                               "CCCCCCCCCCCCCCCCCCCC"_dna5,
                               "TTTTTTTTTTTTTTTTTTTT"_dna5,
                               "GGGGGGGGGGGGGGGGGGGG"_dna5 };
    

    std::string base_prefix{"classify_base_build"};
    GanonBuild::Config cfg_build;
    cfg_build.bin_size_bits      = 5000;
    cfg_build.quiet              = true;
    cfg_build.kmer_size          = 10;
    cfg_build.output_filter_file = base_prefix + ".ibf";
    cfg_build.reference_files = aux::write_sequences_files( base_prefix, "fasta", seqs, ids );
    REQUIRE( GanonBuild::run( cfg_build ) );

    SECTION( "with default config." )
    {
        std::string prefix{"classify_base_build"};
        auto cfg = config_classify::defaultConfig( prefix );
        cfg.ibf ={ base_prefix + ".ibf" };
        aux::write_sequences( "Areads.fasta", { "AAAAAAAAAAAAAA"_dna5 }, { "readA" });
        cfg.single_reads ={ "Areads.fasta" };
        
        REQUIRE( GanonClassify::run( cfg ) );
        config_classify::Res res{cfg};
        config_classify::sanity_check(cfg, res);

        REQUIRE( res.all["readA"]["0"] == 5 );
        REQUIRE( res.all["readA"]["1"] == 0 );
        REQUIRE( res.all["readA"]["2"] == 5 );
        REQUIRE( res.all["readA"]["3"] == 0 );
        
    }
}

SCENARIO( "classifying reads with errors", "[ganon-classify]" )
{
    SECTION( "..." )
    {
        REQUIRE( true);
    }
}