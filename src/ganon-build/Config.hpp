#pragma once

#include <cinttypes>
#include <string>
#include <vector>

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
    bool        update_complete;

    // General options
    std::vector< std::string > reference_files;
    std::uint16_t              threads;
    std::uint16_t              build_threads;
    bool                       verbose;
};
