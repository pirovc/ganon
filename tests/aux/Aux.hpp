#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <cereal/archives/binary.hpp>

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

inline seqan3::interleaved_bloom_filter<> load_ibf( const std::string& file ){
    seqan3::interleaved_bloom_filter<> filter;
    std::ifstream              is( file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( filter );
    return filter;
}

} // namespace aux
