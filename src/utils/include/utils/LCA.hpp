#pragma once

#include <robin_hood.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>

class LCA
{
public:
    LCA() = default;

    void        addEdge( const std::string& father, const std::string& son );
    void        doEulerWalk( const std::string& root_node );
    int         getLCA( unsigned int u, unsigned int v ) const;
    std::string getLCA( const std::vector< std::string >& taxIds ) const;

private:
    void depthFirstSearch( const std::string& current, unsigned int depth );
    void preProcessRMQ();
    int  queryRMQ( unsigned int i, unsigned int j ) const;

    static constexpr int first_appearance_init = -1;

    robin_hood::unordered_map< std::string, std::vector< std::string > > m_parents;
    std::vector< int >                                                   m_euler;
    std::vector< int >                                                   m_depth;
    std::vector< int >                                                   m_firstAppearance;
    int                                                                  m_vertices = 0;
    robin_hood::unordered_map< std::string, unsigned int >               m_encode;
    std::vector< std::string >                                           m_decode;
    std::vector< std::vector< int > >                                    m_M;
};

inline void LCA::addEdge( const std::string& father, const std::string& son )
{
    if ( m_encode.count( father ) == 0 )
    {
        m_encode.insert( { father, m_vertices } );
        m_decode.insert( m_decode.begin() + m_vertices, father );
        ++m_vertices;
    }

    if ( m_encode.count( son ) == 0 )
    {
        m_encode.insert( { son, m_vertices } );
        m_decode.insert( m_decode.begin() + m_vertices, son );
        ++m_vertices;
    }

    if ( m_parents.count( father ) == 0 )
    {
        m_parents[father] = { son };
    }
    else
    {
        m_parents[father].emplace_back( son );
    }
}

inline void LCA::depthFirstSearch( const std::string& current, unsigned int depth )
{
    const auto currentEncoded = m_encode[current];

    // marking first appearance for current node
    if ( m_firstAppearance[currentEncoded] == first_appearance_init )
    {
        m_firstAppearance[currentEncoded] = m_euler.size();
    }

    // pushing root to euler walk
    m_euler.push_back( currentEncoded );
    // pushing depth of current node
    m_depth.push_back( depth );

    for ( const auto& node : m_parents[current] )
    {
        depthFirstSearch( node, depth + 1 );
        m_euler.push_back( currentEncoded );
        m_depth.push_back( depth );
    }
}

inline void LCA::doEulerWalk( const std::string& root_node )
{
    m_firstAppearance.resize( m_vertices, first_appearance_init );
    depthFirstSearch( root_node, 0 );
    preProcessRMQ();
}

// <O(N logN) Preprocessing time, O(1) Query time>
inline void LCA::preProcessRMQ()
{
    const auto size     = m_depth.size();
    const int  logDepth = std::ceil( std::log2( m_depth.size() ) );

    m_M.resize( size, std::vector< int >( logDepth ) );

    for ( auto i = 0u; i < size; ++i )
    {
        m_M[i].front() = i;
    }

    // compute values from smaller to bigger intervals
    for ( unsigned int j = 1; 1u << j <= size; j++ )
    {
        for ( unsigned int i = 0; i + ( 1 << j ) - 1 < size; i++ )
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

inline int LCA::queryRMQ( unsigned int i, unsigned int j ) const
{
    if ( i > j )
    {
        std::swap( i, j );
    }

    const auto k     = static_cast< int >( std::log2( j - i + 1 ) );
    const auto term1 = m_M[i][k];
    const auto term2 = m_M[j - ( 1 << k ) + 1][k];

    if ( m_depth[term1] <= m_depth[term2] )
    {
        return term1;
    }
    else
    {
        return term2;
    }
}

inline int LCA::getLCA( unsigned int u, unsigned int v ) const
{
    assert( u < m_firstAppearance.size() );
    assert( v < m_firstAppearance.size() );

    // trivial case
    if ( u == v )
    {
        return u;
    }

    if ( m_firstAppearance[u] > m_firstAppearance[v] )
    {
        std::swap( u, v );
    }

    // doing RMQ in the required range
    return m_euler[queryRMQ( m_firstAppearance[u], m_firstAppearance[v] )];
}

inline std::string LCA::getLCA( const std::vector< std::string >& taxIds ) const
{
    assert( taxIds.size() > 1 ); // Ideally should return itself if size==1

    int lca = getLCA( m_encode.at( taxIds[0] ), m_encode.at( taxIds[1] ) );
    for ( unsigned int i = 2; i < taxIds.size(); ++i )
        lca = getLCA( lca, m_encode.at( taxIds[i] ) );

    return m_decode.at( lca );
}