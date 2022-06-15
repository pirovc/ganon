#pragma once

#include <cereal/archives/binary.hpp>

struct IBFConfig
{
    IBFConfig()
    : n_bins{ 0 }
    , max_hashes_bin{ 0 }
    , hash_functions{ 0 }
    , kmer_size{ 0 }
    , window_size{ 0 }
    , bin_size_bits{ 0 }
    , max_fp{ 0 }
    {
    }

    uint64_t n_bins;
    uint64_t max_hashes_bin;
    uint8_t  hash_functions;
    uint8_t  kmer_size;
    uint16_t window_size;
    uint64_t bin_size_bits;
    double   max_fp;
    double   true_max_fp;
    double   true_avg_fp;

    template < class Archive >
    void serialize( Archive& archive )
    {
        archive( n_bins,
                 max_hashes_bin,
                 hash_functions,
                 kmer_size,
                 window_size,
                 bin_size_bits,
                 max_fp,
                 true_max_fp,
                 true_avg_fp );
    }
};

inline std::ostream& operator<<( std::ostream& stream, const IBFConfig& ibf_config )
{
    constexpr auto newl{ "\n" };
    constexpr auto separator{ "----------------------------------------------------------------------" };

    stream << "ibf_config:" << newl;
    stream << "n_bins         " << ibf_config.n_bins << newl;
    stream << "max_hashes_bin " << ibf_config.max_hashes_bin << newl;
    stream << "hash_functions " << unsigned( ibf_config.hash_functions ) << newl;
    stream << "kmer_size      " << unsigned( ibf_config.kmer_size ) << newl;
    stream << "window_size    " << ibf_config.window_size << newl;
    stream << "bin_size_bits  " << ibf_config.bin_size_bits << newl;
    stream << "max_fp         " << ibf_config.max_fp << newl;
    stream << "true_max_fp    " << ibf_config.true_max_fp << newl;
    stream << "true_avg_fp    " << ibf_config.true_avg_fp << newl;
    stream << separator << newl;
    return stream;
}
