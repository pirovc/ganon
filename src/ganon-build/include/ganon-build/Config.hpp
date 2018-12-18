#pragma once

#include <cinttypes>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

namespace GanonBuild
{

struct Config
{
    static constexpr std::uint64_t MBinBits = 8388608;

    // Build options
    std::string   seqid_bin_file;
    std::string   output_filter_file;
    std::uint64_t filter_size;
    std::uint16_t kmer_size;
    std::uint16_t hash_functions;

    // Update options
    std::string update_filter_file;
    bool        update_complete = false;

    // General options
    std::vector< std::string > reference_files;
    std::uint16_t              threads;
    std::uint16_t              build_threads;
    bool                       verbose;
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };

    stream << newl;
    stream << "--seqid-bin-file      " << config.seqid_bin_file << newl;
    stream << "--output-filter-file  " << config.output_filter_file << newl;
    stream << "--update-filter-file  " << config.update_filter_file << newl;
    stream << "--reference-files     " << newl;

    for ( const auto& file : config.reference_files )
    {
        stream << "                      " << file << newl;
    }

    stream << "--update-complete     " << config.update_complete << newl;

    stream << "--filter-size         " << std::fixed << std::setprecision( 2 )
           << static_cast< float >( config.filter_size ) / static_cast< float >( Config::MBinBits ) << newl;

    stream << "--filter-size-bits    " << config.filter_size << newl;
    stream << "--kmer-size           " << config.kmer_size << newl;
    stream << "--hash-functions      " << config.hash_functions << newl;
    stream << "--threads             " << config.threads << newl;
    stream << "--verbose             " << config.verbose << newl;
    stream << std::endl;

    return stream;
}

} // namespace GanonBuild
