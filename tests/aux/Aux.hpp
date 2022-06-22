#pragma once

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <utils/IBFConfig.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>


namespace aux
{

inline int fileLines( const std::string& file )
{
    int           count = 0;
    std::ifstream f( file );
    std::string   line;
    while ( getline( f, line ) )
        count++;
    return count;
}

inline unsigned int fileSizeBytes( const std::string& file )
{
    std::ifstream f( file, std::ios::binary | std::ios::ate );
    return f.tellg();
}

inline bool filesAreEqual( const std::string& file1, const std::string& file2 )
{
    std::ifstream     stream1{ file1 };
    const std::string data1{ std::istreambuf_iterator< char >{ stream1 }, std::istreambuf_iterator< char >{} };

    std::ifstream     stream2{ file2 };
    const std::string data2{ std::istreambuf_iterator< char >{ stream2 }, std::istreambuf_iterator< char >{} };

    return data1 == data2;
}

inline bool filesAreEqualSorted( const std::string& file1, const std::string& file2 )
{
    std::ifstream stream1{ file1 };
    std::string   data1{ std::istreambuf_iterator< char >{ stream1 }, std::istreambuf_iterator< char >{} };

    std::ifstream stream2{ file2 };
    std::string   data2{ std::istreambuf_iterator< char >{ stream2 }, std::istreambuf_iterator< char >{} };

    std::sort( data1.begin(), data1.end() );
    std::sort( data2.begin(), data2.end() );

    return data1 == data2;
}

inline bool fileIsEmpty( const std::string& file )
{
    std::ifstream stream{ file };
    if ( stream.peek() == std::ifstream::traits_type::eof() )
        return true;
    else
        return false;
}

inline std::vector< std::vector< std::string > > parse_tsv( const std::string& file )
{
    std::string                               line;
    std::ifstream                             infile( file );
    std::vector< std::vector< std::string > > parsed;

    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field );
        parsed.push_back( fields );
    }
    infile.close();
    return parsed;
}

inline seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed > load_ibf(
    const std::string& file, IBFConfig& ibf_config, std::vector< std::tuple< uint64_t, std::string > >& bin_map )
{
    std::ifstream              is( file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );

    std::tuple< int, int, int >                                           parsed_version;
    std::vector< std::tuple< std::string, uint64_t > >                    hashes_count_std;
    seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed > filter;

    archive( parsed_version );
    archive( ibf_config );
    archive( hashes_count_std );
    archive( bin_map );
    archive( filter );

    return filter;
}

template < typename T >
inline std::vector< T > vconcat( std::vector< T > v1, std::vector< T > v2 )
{
    std::vector< T > cv{ v1 };
    cv.insert( cv.end(), v2.begin(), v2.end() );
    return cv;
}

} // namespace aux
