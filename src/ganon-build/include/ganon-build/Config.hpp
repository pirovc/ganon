#pragma once

#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include <stdio.h>
#include <dirent.h>

namespace GanonBuild
{

struct Config
{

public:
    static constexpr uint32_t MBinBits = 8388608;

    // Required
    std::string                seqid_bin_file     = "";
    std::string                output_filter_file = "";
    std::vector< std::string > reference_files;

    // Default
    uint32_t    filter_size        = 0;
    uint64_t    filter_size_bits   = 0;
    uint16_t    kmer_size          = 19;
    uint16_t    hash_functions     = 3;
    std::string update_filter_file = "";
    bool        update_complete    = false;
    uint16_t    threads            = 2;
    bool        verbose            = false;
    bool        quiet              = false;
    std::string directory_reference_files = "";
    std::string extension = "";

    // hidden
    uint32_t n_batches = 1000;
    uint32_t n_refs    = 400;

    // private
    const uint16_t min_threads = 2;
    uint16_t       build_threads;

    bool hasEnding (std::string const &fullString, std::string const &ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
        } else {
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
        if (!directory_reference_files.empty() && !extension.empty())
        {
            struct dirent *entry = nullptr;
            DIR *dp = nullptr;
            dp = opendir(directory_reference_files.c_str());
            if (dp != nullptr) {
                while ((entry = readdir(dp))){
                    if (hasEnding(entry->d_name, extension))
                        reference_files.push_back(directory_reference_files + "/" + entry->d_name);
                }
            }
            closedir(dp);
        }

        if (reference_files.empty()){
            std::cerr << "Please provide reference sequence files with the parameters --reference-files or/and with --directory-reference-files and --extension" << std::endl;
            return false;
        }

        if ( threads < min_threads )
            threads = min_threads;

        build_threads = threads - 1; // -1 reading files

        // Skip variables if updating, loads from existing filter file
        if ( !update_filter_file.empty() )
        {
            std::cerr << "WARNING: --filter-size[-bits], --kmer-size --hash-funtions ignored, using metadata from "
                         "--update-filter-file"
                      << std::endl;
            kmer_size        = 0;
            hash_functions   = 0;
            filter_size      = 0;
            filter_size_bits = 0;
        }
        else
        {
            if ( filter_size_bits == 0 )
            {
                if ( filter_size == 0 )
                {
                    std::cerr << "--filter-size or --filter-size-bits are required" << std::endl;
                    return false;
                }
                else
                {
                    filter_size_bits = filter_size * MBinBits;
                }
            }
            else
            {
                filter_size = filter_size_bits / MBinBits;
            }
        }

        if ( n_batches < 1 )
            n_batches = 1;

        if ( n_refs < 1 )
            n_refs = 1;

        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };

    stream << newl;
    stream << "--seqid-bin-file      " << config.seqid_bin_file << newl;
    stream << "--output-filter-file  " << config.output_filter_file << newl;
    stream << "--filter-size         " << config.filter_size << newl;
    stream << "--filter-size-bits    " << config.filter_size_bits << newl;
    stream << "--hash-functions      " << config.hash_functions << newl;
    stream << "--kmer-size           " << config.kmer_size << newl;
    stream << "--n-batches           " << config.n_batches << newl;
    stream << "--n-refs              " << config.n_refs << newl;
    stream << "--threads             " << config.threads << newl;
    stream << "--update-filter-file  " << config.update_filter_file << newl;
    stream << "--update-complete     " << config.update_complete << newl;
    stream << "--verbose             " << config.verbose << newl;
    stream << "--reference-files     " << newl;
    for ( const auto& file : config.reference_files )
    {
        stream << "                      " << file << newl;
    }
    stream << std::endl;

    return stream;
}

} // namespace GanonBuild
