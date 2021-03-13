#pragma once

#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include <dirent.h>
#include <stdio.h>

namespace GanonBuild
{

struct Config
{

public:
    std::vector< std::string > reference_files;
    std::string                directory_reference_files = "";
    std::string                extension                 = "";

    std::string seqid_bin_file     = "";
    std::string output_filter_file = "";
    std::string update_filter_file = "";
    bool        update_complete    = false;

    uint32_t filter_size_mb = 0;
    uint64_t bin_size_bits  = 0;

    uint8_t  kmer_size      = 19;
    uint16_t hash_functions = 3;

    uint16_t threads   = 2;
    uint32_t n_refs    = 400;
    uint32_t n_batches = 1000;
    bool     verbose   = false;
    bool     quiet     = false;

    uint16_t threads_build;

    bool hasEnding( std::string const& fullString, std::string const& ending )
    {
        if ( fullString.length() >= ending.length() )
        {
            return ( 0 == fullString.compare( fullString.length() - ending.length(), ending.length(), ending ) );
        }
        else
        {
            return false;
        }
    }

    bool validate()
    {

        if ( seqid_bin_file.empty() || output_filter_file.empty() )
        {
            std::cerr << "--seqid-bin-file and --output-filter-file are mandatory" << std::endl;
            return false;
        }

        // add references from folder
        if ( !directory_reference_files.empty() && !extension.empty() )
        {
            struct dirent* entry = nullptr;
            DIR*           dp    = nullptr;
            dp                   = opendir( directory_reference_files.c_str() );
            if ( dp != nullptr )
            {
                while ( ( entry = readdir( dp ) ) )
                {
                    if ( hasEnding( entry->d_name, extension ) )
                        reference_files.push_back( directory_reference_files + "/" + entry->d_name );
                }
            }
            closedir( dp );
        }

        if ( reference_files.empty() )
        {
            std::cerr << "Please provide reference sequence files with the parameters --reference-files or/and with "
                         "--directory-reference-files and --extension"
                      << std::endl;
            return false;
        }

        if ( threads <= 2 )
        {
            threads_build = 1;
        }
        else
        {
            threads_build = threads - 1; // extra reading thread
        }

        if ( n_batches < 1 )
            n_batches = 1;

        if ( n_refs < 1 )
            n_refs = 1;

        // Skip variables if updating, loads from existing filter file
        if ( !update_filter_file.empty() )
        {
            if ( verbose && ( filter_size_mb > 0 || bin_size_bits > 0 ) )
            {
                std::cerr << "WARNING: --filter-size-mb and --bin-size-bits ignored when updating" << std::endl;
            }
            filter_size_mb = 0;
            bin_size_bits  = 0;
        }
        else
        {
            if ( bin_size_bits == 0 && filter_size_mb == 0 )
            {
                std::cerr << "--filter-size-mb or --bin-size-bits are required" << std::endl;
                return false;
            }
        }

        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };
    constexpr auto separator{ "----------------------------------------------------------------------" };

    stream << separator << newl;
    stream << "--reference-files     " << newl;
    for ( const auto& file : config.reference_files )
    {
        stream << "                      " << file << newl;
    }
    stream << "--seqid-bin-file      " << config.seqid_bin_file << newl;
    stream << "--output-filter-file  " << config.output_filter_file << newl;
    stream << "--update-filter-file  " << config.update_filter_file << newl;
    stream << "--update-complete     " << config.update_complete << newl;
    if ( config.bin_size_bits > 0 )
        stream << "--bin-size-bits       " << config.bin_size_bits << newl;
    if ( config.filter_size_mb > 0 )
        stream << "--filter-size-mb      " << config.filter_size_mb << newl;
    stream << "--hash-functions      " << config.hash_functions << newl;
    stream << "--kmer-size           " << config.kmer_size << newl;
    stream << "--threads             " << config.threads << newl;
    stream << "--n-refs              " << config.n_refs << newl;
    stream << "--n-batches           " << config.n_batches << newl;
    stream << "--verbose             " << config.verbose << newl;
    stream << "--quiet               " << config.quiet << newl;
    stream << separator << newl;

    return stream;
}

} // namespace GanonBuild
