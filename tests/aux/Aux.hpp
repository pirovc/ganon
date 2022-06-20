#pragma once

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <utils/IBFConfig.hpp>
#include <utils/dna4_traits.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

using sequences_type       = std::vector< seqan3::dna4_vector >;
using sequence_record_type = seqan3::sequence_record< seqan3::type_list< std::vector< seqan3::dna4 >, std::string >,
                                                      seqan3::fields< seqan3::field::seq, seqan3::field::id > >;

namespace aux
{

inline void write_sequences( const std::string file, const sequences_type& seqs, const std::vector< std::string >& ids )
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
                                                         const sequences_type& seqs )
{
    std::vector< std::string > output_files;
    for ( uint16_t i = 0; i < seqs.size(); ++i )
    {
        auto        id = std::to_string( i );
        std::string filename{ prefix + id + "." + suffix };
        write_sequences( filename, { seqs[i] }, { id } );
        output_files.push_back( filename );
    }
    return output_files;
}

inline std::vector< std::string > write_sequences_files( const std::string                 prefix,
                                                         const std::string                 suffix,
                                                         const sequences_type&             seqs,
                                                         const std::vector< std::string >& ids )
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

inline void write_input_file( std::string out_file, std::vector< std::string >& files )
{
    std::ofstream output_file{ out_file };
    for ( auto& file : files )
    {
        output_file << std::filesystem::canonical( file ).c_str() << '\n';
    }
    output_file.close();
}

inline void write_input_file_target( std::string                 out_file,
                                     std::vector< std::string >& files,
                                     std::vector< std::string >& targets_file )
{
    std::ofstream output_file{ out_file };
    size_t        i = 0;
    for ( auto& file : files )
    {
        output_file << std::filesystem::canonical( file ).c_str() << '\t' << targets_file[i] << '\n';
        i++;
    }
    output_file.close();
}

inline void write_input_file_seqs( std::string                 out_file,
                                   std::vector< std::string >& files,
                                   std::vector< std::string >& targets_seq )
{
    std::ofstream output_file{ out_file };
    size_t        i = 0;
    for ( auto& file : files )
    {
        seqan3::sequence_file_input< raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq > > fin{
            file
        };
        for ( auto const& [header, seq] : fin )
        {
            output_file << std::filesystem::canonical( file ).c_str() << '\t' << targets_seq[i] << '\t' << header
                        << '\n';
            i++;
        }
    }
    output_file.close();
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
