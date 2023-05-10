#include "aux/Aux.hpp"

#include <utils/LCA.hpp>

#include <catch2/catch.hpp>

LCA pre_process_lca( std::string file )
{
    LCA lca;
    for ( auto const& line : aux::parse_tsv( file ) )
        lca.addEdge( line[1], line[0] );
    lca.doEulerWalk( "1" );
    return lca;
}

// lca_values[lca] = {val1, val2, ..., valn}
typedef std::unordered_map< std::string, std::vector< std::string > > TLcaValues;

SCENARIO( "LCA general test", "[utils][lca]" )
{
    GIVEN( "A pre-processed LCA" )
    {
        LCA lca = pre_process_lca( "lca/tree.tax" );
        WHEN( "LCA values are requested" )
        {
            TLcaValues lca_values;
            lca_values["D0"] = { "E0", "E1" };
            lca_values["C3"] = { "C3", "F4" };
            lca_values["A0"] = { "G0", "C3", "D5" };
            lca_values["1"]  = { "G0", "G5" };

            THEN( "they are valid" )
            {
                for ( auto& v : lca_values )
                    REQUIRE( lca.getLCA( v.second ) == v.first );
            }
        }
    }
}

// SCENARIO( "Invalid nodes", "[utils][lca]" )
// {
//     GIVEN( "A pre-processed LCA" )
//     {
//         LCA lca = pre_process_lca( "lca/tree.tax" );
//         WHEN( "LCA values are requested" )
//         {
//             TLcaValues lca_values;
//             lca_values["1"] = { "xxxxx", "aaaaaa" };
//             lca_values["1"] = { "xxxxx", "pppppp", "xccccc" };
//             lca_values["1"] = { "xxxxx", "pppppp", "xccccc", "E2" };
//             THEN( "they are valid" )
//             {
//                 for ( auto& v : lca_values )
//                     REQUIRE( lca.getLCA( v.second ) == v.first );
//             }
//         }
//     }
// }

SCENARIO( "Reverse test", "[utils][lca]" )
{
    GIVEN( "A pre-processed LCA" )
    {
        LCA lca = pre_process_lca( "lca/tree.tax" );

        WHEN( "LCA values are requested" )
        {
            TLcaValues lca_values;
            lca_values["B1"] = { "B1", "C2" };
            lca_values["B1"] = { "C2", "B1" };
            lca_values["B0"] = { "C0", "E1", "F2" };
            lca_values["B0"] = { "C0", "F2", "E1" };
            lca_values["B0"] = { "F2", "C0", "E1" };
            lca_values["B0"] = { "F2", "E1", "C0" };
            lca_values["B0"] = { "E1", "F2", "C0" };
            lca_values["B0"] = { "E1", "C0", "F2" };

            THEN( "they are valid" )
            {
                for ( auto& v : lca_values )
                    REQUIRE( lca.getLCA( v.second ) == v.first );
            }
        }
    }
}

SCENARIO( "NCBI test", "[utils][lca]" )
{
    GIVEN( "A pre-processed LCA" )
    {
        LCA lca = pre_process_lca( "lca/ncbi.tax" );
        WHEN( "LCA values are requested" )
        {
            TLcaValues lca_values;
            // Bacteria
            lca_values["1224"] = { "366602", "470" };
            lca_values["2"]    = { "366602", "470", "1406" };
            // Euryarchaeota
            lca_values["2290931"] = { "2223", "51589" };
            // Virus
            lca_values["10239"] = { "2025595", "491893" };

            // mix
            lca_values["1"] = { "366602", "2223" };
            lca_values["1"] = { "51589", "2025595" };
            lca_values["1"] = { "366602", "470", "1406", "2223", "51589", "2025595", "491893" };

            THEN( "they are valid" )
            {
                for ( auto& v : lca_values )
                    REQUIRE( lca.getLCA( v.second ) == v.first );
            }
        }
    }
}