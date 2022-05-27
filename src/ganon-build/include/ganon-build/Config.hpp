#pragma once

#include <filesystem>
#include <iostream>

namespace GanonBuild
{

struct Config
{

public:
    std::string input_file;
    std::string output_file;
    std::string tmp_output_folder = "";
    double      max_fp            = 0.05;
    double      filter_size       = 0;
    uint8_t     kmer_size         = 19;
    uint16_t    window_size       = 32;
    uint8_t     hash_functions    = 0;
    uint16_t    threads           = 2;
    uint32_t    n_refs            = 400;
    uint32_t    n_batches         = 1000;
    bool        verbose           = false;
    bool        quiet             = false;

    uint16_t threads_build = 1;

    bool validate()
    {

        if ( input_file.empty() )
        {
            if ( !quiet )
                std::cerr << "--input-file is mandatory" << std::endl;
            return false;
        }
        else
        {
            if ( !std::filesystem::exists( input_file ) )
            {
                if ( !quiet )
                    std::cerr << "file not found: " << input_file << std::endl;
                return false;
            }
            else if ( std::filesystem::file_size( input_file ) == 0 )
            {
                if ( !quiet )
                    std::cerr << "file is empty: " << input_file << std::endl;
                return false;
            }
        }

        if ( output_file.empty() )
        {
            if ( !quiet )
                std::cerr << "--output-file is mandatory" << std::endl;
            return false;
        }

        // 2 threads default: one read input, one process
        if ( threads <= 2 )
            threads_build = 1;
        else
            threads_build = threads - 1;

        if ( n_batches < 1 )
            n_batches = 1;

        if ( n_refs < 1 )
            n_refs = 1;

        // Skip variables if updating, loads from existing filter file
        if ( filter_size == 0 && max_fp == 0 )
        {
            if ( !quiet )
                std::cerr << "--false-positive or --filter-size are mandatory" << std::endl;
        }

        if ( window_size > 0 && window_size < kmer_size )
        {
            if ( !quiet )
                std::cerr << "--window-size has to be >= --kmer-size" << std::endl;
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
    stream << "--input-file        " << config.input_file << newl;
    stream << "--output-file       " << config.output_file << newl;
    stream << "--tmp-output-folder " << config.tmp_output_folder << newl;
    stream << "--max-fp            " << config.max_fp << newl;
    stream << "--filter-size       " << config.filter_size << newl;
    stream << "--kmer-size         " << unsigned( config.kmer_size ) << newl;
    stream << "--window-size       " << config.window_size << newl;
    stream << "--hash-functions    " << unsigned( config.hash_functions ) << newl;
    stream << "--threads           " << config.threads << newl;
    stream << "--n-refs            " << config.n_refs << newl;
    stream << "--n-batches         " << config.n_batches << newl;
    stream << "--verbose           " << config.verbose << newl;
    stream << "--quiet             " << config.quiet << newl;
    stream << separator << newl;

    return stream;
}

} // namespace GanonBuild
