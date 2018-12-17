#pragma once

#include <cinttypes>
#include <iomanip>
#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace GanonClassify
{

struct Config
{
    std::string                                                                    output_file;
    std::string                                                                    output_unclassified_file;
    std::string                                                                    filter_hierarchy;
    uint16_t                                                                       max_error;
    uint16_t                                                                       threads = 3;
    uint16_t                                                                       clas_threads;
    bool                                                                           output_unclassified;
    bool                                                                           unique_filtering;
    int16_t                                                                        max_error_unique = -1;
    std::vector< std::string >                                                     bloom_filter_files;
    std::vector< std::string >                                                     group_bin_files;
    std::vector< std::string >                                                     reads;
    std::map< std::string, std::vector< std::tuple< std::string, std::string > > > filters;
    bool                                                                           verbose;
    bool                                                                           testing = false; // internal
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };

    stream << newl;
    stream << "--filter-hierarchy          " << config.filter_hierarchy << newl;
    stream << "--bloom-filter,--group-bin  " << newl;
    for ( auto const& hierarchy : config.filters )
    {
        for ( auto const& file : hierarchy.second )
        {
            stream << "                            (" << hierarchy.first << ") " << std::get< 0 >( file ) << ", "
                   << std::get< 1 >( file ) << newl;
        }
    }
    stream << "--max-error                 " << config.max_error << newl;
    stream << "--max-error-unique          " << config.max_error_unique << newl;
    stream << "--output-file               " << config.output_file << newl;
    stream << "--output-unclassified-file  " << config.output_unclassified_file << newl;
    stream << "--verbose                   " << config.verbose << newl;
    stream << "--threads                   " << config.threads << newl;
    stream << "--reads                     " << newl;
    for ( const auto& s : config.reads )
        stream << "                            " << s << newl;
    stream << newl;

    return stream;
}

} // namespace GanonClassify
