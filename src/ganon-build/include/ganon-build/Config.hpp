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
    std::string mode              = "avg";
    double      max_fp            = 0.05;
    double      filter_size       = 0;
    uint8_t     kmer_size         = 19;
    uint16_t    window_size       = 31;
    uint8_t     hash_functions    = 0;
    uint64_t    min_length        = 0;
    uint16_t    threads           = 1;
    bool        verbose           = false;
    bool        quiet             = false;

    const uint8_t max_hash_functions = 5;

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
                    std::cerr << "--input-file not found: " << input_file << std::endl;
                return false;
            }
            else if ( std::filesystem::file_size( input_file ) == 0 )
            {
                if ( !quiet )
                    std::cerr << "--input-file is empty: " << input_file << std::endl;
                return false;
            }
        }

        if ( output_file.empty() )
        {
            if ( !quiet )
                std::cerr << "--output-file is mandatory" << std::endl;
            return false;
        }

        if ( tmp_output_folder != "" && !std::filesystem::exists( tmp_output_folder ) )
        {
            if ( !quiet )
                std::cerr << "--tmp-output-folder not found" << std::endl;
            return false;
        }

        if ( hash_functions > max_hash_functions )
        {
            if ( !quiet )
                std::cerr << "--hash-functions must be <=5" << std::endl;
            return false;
        }

        if ( filter_size == 0 && max_fp == 0 )
        {
            if ( !quiet )
                std::cerr << "--max-fp or --filter-size is mandatory" << std::endl;
            return false;
        }

        if ( filter_size > 0 )
        {
            max_fp = 0;
        }

        if ( window_size < kmer_size )
        {
            if ( !quiet )
                std::cerr << "--window-size has to be >= --kmer-size" << std::endl;
            return false;
        }

        if ( mode != "avg" && mode != "smaller" && mode != "smallest" && mode != "faster" && mode != "fastest" )
        {
            if ( !quiet )
                std::cerr << "Invalid --mode" << std::endl;
            return false;
        }

        if ( kmer_size > 32 )
        {
            if ( !quiet )
                std::cerr << "--kmer-size has to be <= 32" << std::endl;
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
    stream << "--mode              " << config.mode << newl;
    stream << "--min-length        " << config.min_length << newl;
    stream << "--threads           " << config.threads << newl;
    stream << "--verbose           " << config.verbose << newl;
    stream << "--quiet             " << config.quiet << newl;
    stream << separator << newl;

    return stream;
}

} // namespace GanonBuild
