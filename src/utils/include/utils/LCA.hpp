#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// required: root node == "1" and father == "0"
class LCA
{
private:
    unordered_map< string, vector< string > > parents;
    vector< int >                             euler;
    vector< int >                             depth;
    unique_ptr< int[] >                       firstAppearance;
    int                                       vertices;
    unordered_map< string, int >              encode;
    unordered_map< int, string >              decode;
    unique_ptr< unique_ptr< int[] >[] >       M;
    void                                      depthFirstSearch( string current, int depth );
    void                                      preProcessRMQ();
    int                                       queryRMQ( int i, int j );

public:
    LCA();
    void   addEdge( string father, string son );
    void   doEulerWalk();
    int    getLCA( int u, int v );
    string getLCA( vector< string >& taxIds );
};

inline LCA::LCA()
{
    vertices = 0;
}

inline void LCA::addEdge( string father, string son )
{
    if ( encode.count( father ) == 0 )
    {
        encode.insert( { father, vertices } );
        decode.insert( { vertices, father } );
        vertices++;
    }
    if ( encode.count( son ) == 0 )
    {
        encode.insert( { son, vertices } );
        decode.insert( { vertices, son } );
        vertices++;
    }
    if ( parents.count( father ) == 0 )
    {
        vector< string > children;
        children.push_back( son );
        parents[father] = children;
    }
    else
    {
        parents[father].push_back( son );
    }
}

inline void LCA::depthFirstSearch( string current, int _depth )
{
    // marking first appearance for current node
    if ( firstAppearance[encode[current]] == -1 )
    {
        firstAppearance[encode[current]] = euler.size();
    }
    // pushing root to euler walk
    euler.push_back( encode[current] );
    // pushing depth of current node
    this->depth.push_back( _depth );
    for ( unsigned int i = 0; i < parents[current].size(); i++ )
    {
        depthFirstSearch( parents[current][i], _depth + 1 );
        euler.push_back( encode[current] );
        this->depth.push_back( _depth );
    }
}

inline void LCA::doEulerWalk()
{
    firstAppearance = make_unique< int[] >( vertices );
    for ( int i = 0; i < vertices; i++ )
    {
        firstAppearance[i] = -1;
    }
    depthFirstSearch( "1", 0 );
    preProcessRMQ();
}

// <O(N logN) Preprocessing time, O(1) Query time>
inline void LCA::preProcessRMQ()
{

    M = make_unique< unique_ptr< int[] >[] >( depth.size() );

    int logDepth = ceil( log2( depth.size() ) );
    for ( unsigned int i = 0; i < depth.size(); i++ )
    {
        M[i]    = make_unique< int[] >( logDepth );
        M[i][0] = i; // initialize M for the intervals with length 1
    }


    // compute values from smaller to bigger intervals
    for ( unsigned int j = 1; 1u << j <= depth.size(); j++ )
    {
        for ( unsigned int i = 0; i + ( 1 << j ) - 1 < depth.size(); i++ )
        {
            if ( depth[M[i][j - 1]] < depth[M[i + ( 1 << ( j - 1 ) )][j - 1]] )
            {
                M[i][j] = M[i][j - 1];
            }
            else
            {
                M[i][j] = M[i + ( 1 << ( j - 1 ) )][j - 1];
            }
        }
    }
}

inline int LCA::queryRMQ( int i, int j )
{
    if ( i > j )
    {
        swap( i, j );
    }

    int k = log2( j - i + 1 );

    if ( depth[M[i][k]] <= depth[M[j - ( 1 << k ) + 1][k]] )
    {
        return M[i][k];
    }
    else
    {
        return M[j - ( 1 << k ) + 1][k];
    }
}

inline int LCA::getLCA( int u, int v )
{
    // trivial case
    if ( u == v )
    {
        return u;
    }

    // check for invalid nodes
    if ( u == 0 || v == 0 )
    {
        return 0;
    }

    if ( firstAppearance[u] > firstAppearance[v] )
    {
        swap( u, v );
    }

    // doing RMQ in the required range
    return euler[queryRMQ( firstAppearance[u], firstAppearance[v] )];
}

inline string LCA::getLCA( vector< string >& taxIds )
{
    int lca;
    lca = getLCA( encode[taxIds[0]], encode[taxIds[1]] );
    for ( unsigned int i = 2; i < taxIds.size(); i++ )
    {
        lca = getLCA( lca, encode[taxIds[i]] );
    }
    return decode[lca];
}
