#pragma once

#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

using bins_type            = std::vector< uint16_t >;
using sequences_type       = std::vector< seqan3::dna4_vector >;
using ids_type             = std::vector< std::string >;
using sequence_record_type = seqan3::sequence_record< seqan3::type_list< std::vector< seqan3::dna4 >, std::string >,
                                                      seqan3::fields< seqan3::field::seq, seqan3::field::id > >;

namespace aux
{

inline void write_sequences( const std::string file, const sequences_type& seqs, const ids_type& ids )
{
    seqan3::sequence_file_output fout{ file };
    int                          i = 0;
    for ( auto& seq : seqs )
    {
        sequence_record_type rec{ seq, ids[i] };
        fout.push_back( rec );
        i += 1;
    }
}


inline std::vector< std::string > write_sequences_files( const std::string     prefix,
                                                         const std::string     suffix,
                                                         const sequences_type& seqs,
                                                         const ids_type&       ids )
{
    std::vector< std::string > output_files;
    for ( uint16_t i = 0; i < seqs.size(); ++i )
    {
        std::string filename{ prefix + ids[i] + "." + suffix };
        write_sequences( filename, { seqs[i] }, { ids[i] } );
        output_files.push_back( filename );
    }
    return output_files;
}

inline void write_seqid_bin( std::string file, const sequences_type& seqs, const ids_type& ids, const bins_type& bins )
{
    // generate basic seqid_bin -> every sequence in one bin, no fragmentation
    std::ofstream seqid_bin_file{ file };
    uint16_t      i = 0;
    for ( auto& seq : seqs )
    {
        seqid_bin_file << ids[i] << "\t1\t" << std::ranges::size( seq ) << "\t" << bins[i] << '\n';
        i += 1;
    }
    seqid_bin_file.close();
}

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

inline seqan3::interleaved_bloom_filter<> load_ibf( const std::string& file )
{
    seqan3::interleaved_bloom_filter<> filter;
    std::ifstream                      is( file, std::ios::binary );
    cereal::BinaryInputArchive         archive( is );
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
