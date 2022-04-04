#pragma once

#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <seqan3/std/filesystem>
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
    std::string map                = "";
    std::string output_filter_file = "";
    std::string update_filter_file = "";
    bool        update_complete    = false;

    double   false_positive = 0;
    double   filter_size_mb = 0;
    uint64_t bin_size_bits  = 0;

    uint8_t  kmer_size      = 19;
    uint32_t window_size    = 0;
    uint16_t hash_functions = 3;
    bool     count_hashes   = false;

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

    bool check_files( std::vector< std::string > const& files )
    {
        for ( auto const& file : files )
        {
            if ( !std::filesystem::exists( file ) )
            {
                if ( !quiet )
                    std::cerr << "file not found: " << file << std::endl;
                return false;
            }
            else if ( std::filesystem::file_size( file ) == 0 )
            {
                if ( !quiet )
                    std::cerr << "file is empty: " << file << std::endl;
                return false;
            }
        }
        return true;
    }

    bool validate()
    {

        if ( !check_files( reference_files ) )
            return false;

        if ( output_filter_file.empty() )
        {
            if ( !quiet )
                std::cerr << "--output-filter-file is mandatory" << std::endl;
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
            if ( !quiet )
                std::cerr
                    << "Please provide reference sequence files with the parameters --reference-files or/and with "
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
            if ( !check_files( { update_filter_file } ) )
                return false;

            if ( verbose && !quiet && ( filter_size_mb > 0 || bin_size_bits > 0 || false_positive > 0 ) )
            {
                std::cerr << "WARNING: --false-positive, --filter-size-mb and --bin-size-bits ignored when updating"
                          << std::endl;
            }
            filter_size_mb = 0;
            bin_size_bits  = 0;
        }
        else
        {
            if ( bin_size_bits == 0 && filter_size_mb == 0 && false_positive == 0 )
            {
                if ( !quiet )
                    std::cerr << "--false-positive, --filter-size-mb or --bin-size-bits should be provided"
                              << std::endl;
                return false;
            }
        }
        if ( update_complete )
        {
            if ( update_filter_file.empty() || seqid_bin_file.empty() )
            {
                if ( !quiet )
                    std::cerr << "--update-filter-file and --seqid-bin-file are required to perform --update-complete"
                              << std::endl;
                return false;
            }
        }

        if ( window_size > 0 && window_size < kmer_size )
        {
            if ( !quiet )
                std::cerr << "--window-size has to be >= --kmer-size (or 0 to disable windows)" << std::endl;
            return false;
        }

        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };
    constexpr auto separator{ "----------------------------------------------------------------------" };

    stream << separator << newl;
    stream << "--reference-files     " << config.reference_files.size() << " file(s)" << newl;
    stream << "                      " << config.reference_files[0] << newl;
    if ( config.reference_files.size() > 1 )
    {
        stream << "                      ..." << newl;
        stream << "                      " << config.reference_files.back() << newl;
    }
    stream << "--seqid-bin-file      " << config.seqid_bin_file << newl;
    stream << "--map                 " << config.map << newl;
    stream << "--output-filter-file  " << config.output_filter_file << newl;
    stream << "--update-filter-file  " << config.update_filter_file << newl;
    stream << "--update-complete     " << config.update_complete << newl;
    if ( config.false_positive > 0 )
        stream << "--false-positive      " << config.false_positive << newl;
    if ( config.bin_size_bits > 0 )
        stream << "--bin-size-bits       " << config.bin_size_bits << newl;
    if ( config.filter_size_mb > 0 )
        stream << "--filter-size-mb      " << config.filter_size_mb << newl;
    stream << "--kmer-size           " << unsigned( config.kmer_size ) << newl;
    stream << "--window-size         " << config.window_size << newl;
    stream << "--hash-functions      " << config.hash_functions << newl;
    stream << "--count-hashes        " << config.count_hashes << newl;
    stream << "--threads             " << config.threads << newl;
    stream << "--n-refs              " << config.n_refs << newl;
    stream << "--n-batches           " << config.n_batches << newl;
    stream << "--verbose             " << config.verbose << newl;
    stream << "--quiet               " << config.quiet << newl;
    stream << separator << newl;

    return stream;
}

} // namespace GanonBuild
