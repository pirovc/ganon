#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// required: root node == "1" and father == "0"
class LCA
{
public:
    LCA() = default;

    void        addEdge( std::string father, std::string son );
    void        doEulerWalk();
    int         getLCA( int u, int v );
    std::string getLCA( const std::vector< std::string >& taxIds );

private:
    void depthFirstSearch( std::string current, int depth );
    void preProcessRMQ();
    int  queryRMQ( int i, int j );

    std::unordered_map< std::string, std::vector< std::string > > m_parents;
    std::vector< int >                                            m_euler;
    std::vector< int >                                            m_depth;
    std::vector< int >                                            m_firstAppearance;
    int                                                           m_vertices = 0;
    std::unordered_map< std::string, int >                        m_encode;
    std::unordered_map< int, std::string >                        m_decode;
    std::unique_ptr< std::unique_ptr< int[] >[] >                 m_M;
};

inline void LCA::addEdge( std::string father, std::string son )
{
    if ( m_encode.count( father ) == 0 )
    {
        m_encode.insert( { father, m_vertices } );
        m_decode.insert( { m_vertices, father } );
        m_vertices++;
    }
    if ( m_encode.count( son ) == 0 )
    {
        m_encode.insert( { son, m_vertices } );
        m_decode.insert( { m_vertices, son } );
        m_vertices++;
    }
    if ( m_parents.count( father ) == 0 )
    {
        std::vector< std::string > children;
        children.push_back( son );
        m_parents[father] = children;
    }
    else
    {
        m_parents[father].push_back( son );
    }
}

inline void LCA::depthFirstSearch( std::string current, int depth )
{
    // marking first appearance for current node
    if ( m_firstAppearance[m_encode[current]] == -1 )
    {
        m_firstAppearance[m_encode[current]] = m_euler.size();
    }
    // pushing root to euler walk
    m_euler.push_back( m_encode[current] );
    // pushing depth of current node
    m_depth.push_back( depth );
    for ( unsigned int i = 0; i < m_parents[current].size(); i++ )
    {
        depthFirstSearch( m_parents[current][i], depth + 1 );
        m_euler.push_back( m_encode[current] );
        m_depth.push_back( depth );
    }
}

inline void LCA::doEulerWalk()
{
    m_firstAppearance.resize( m_vertices, -1 );
    depthFirstSearch( "1", 0 );
    preProcessRMQ();
}

// <O(N logN) Preprocessing time, O(1) Query time>
inline void LCA::preProcessRMQ()
{

    m_M = std::make_unique< std::unique_ptr< int[] >[] >( m_depth.size() );

    int logDepth = std::ceil( std::log2( m_depth.size() ) );
    for ( unsigned int i = 0; i < m_depth.size(); i++ )
    {
        m_M[i]    = std::make_unique< int[] >( logDepth );
        m_M[i][0] = i; // initialize M for the intervals with length 1
    }


    // compute values from smaller to bigger intervals
    for ( unsigned int j = 1; 1u << j <= m_depth.size(); j++ )
    {
        for ( unsigned int i = 0; i + ( 1 << j ) - 1 < m_depth.size(); i++ )
        {
            if ( m_depth[m_M[i][j - 1]] < m_depth[m_M[i + ( 1 << ( j - 1 ) )][j - 1]] )
            {
                m_M[i][j] = m_M[i][j - 1];
            }
            else
            {
                m_M[i][j] = m_M[i + ( 1 << ( j - 1 ) )][j - 1];
            }
        }
    }
}

inline int LCA::queryRMQ( int i, int j )
{
    if ( i > j )
    {
        std::swap( i, j );
    }

    int k = std::log2( j - i + 1 );

    if ( m_depth[m_M[i][k]] <= m_depth[m_M[j - ( 1 << k ) + 1][k]] )
    {
        return m_M[i][k];
    }
    else
    {
        return m_M[j - ( 1 << k ) + 1][k];
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

    if ( m_firstAppearance[u] > m_firstAppearance[v] )
    {
        std::swap( u, v );
    }

    // doing RMQ in the required range
    return m_euler[queryRMQ( m_firstAppearance[u], m_firstAppearance[v] )];
}

inline std::string LCA::getLCA( const std::vector< std::string >& taxIds )
{
    assert( taxIds.size() > 1 );

    const auto lca = std::accumulate(
        std::next( taxIds.begin(), 2 ),
        taxIds.end(),
        getLCA( m_encode[taxIds[0]], m_encode[taxIds[1]] ),
        [&]( const auto prevLCA, const std::string nextId ) { return getLCA( prevLCA, m_encode[nextId] ); } );

    return m_decode[lca];
}
