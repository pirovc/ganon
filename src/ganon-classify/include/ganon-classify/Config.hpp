#pragma once

#include <algorithm>
#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace GanonClassify
{

struct FilterConfig
{
    std::string bloom_filter_file;
    std::string group_bin_file;
    int16_t     max_error;
    float       min_kmers;
};

struct HierarchyConfig
{
    std::vector< FilterConfig > filters;
    int16_t                     max_error_unique;
    std::string                 output_file;
};

struct Config
{

public:
    // Required
    std::vector< std::string > bloom_filter_files;
    std::vector< std::string > group_bin_files;
    std::vector< std::string > reads_single;
    std::vector< std::string > reads_paired;

    // Defaults
    float                      default_min_kmers = 0.25;
    std::vector< float >       min_kmers;
    std::vector< int16_t >     max_error;
    std::vector< int16_t >     max_error_unique{ -1 };
    std::vector< std::string > filter_hierarchy{ "1" };
    int16_t                    paired_mode                 = -1;
    std::string                output_file                 = "";
    std::string                output_unclassified_file    = "";
    uint16_t                   offset                      = 1;
    uint16_t                   threads                     = 3;
    bool                       verbose                     = false;
    bool                       split_output_file_hierarchy = false;

    // hidden
    uint32_t n_batches = 1000;
    uint32_t n_reads   = 400;

    // private:
    const uint16_t                           min_threads = 3;
    uint16_t                                 clas_threads;
    bool                                     output_unclassified = false;
    std::map< std::string, HierarchyConfig > h_filters;


    bool validate()
    {

        if ( bloom_filter_files.size() == 0 || group_bin_files.size() == 0
             || ( reads_paired.size() == 0 && reads_single.size() == 0 ) )
        {
            std::cerr << "--bloom-filter, --group-bin, --[single|paired]-reads are mandatory" << std::endl;
            return false;
        }
        else if ( reads_paired.size() % 2 != 0 )
        {
            std::cerr << "--paired-reads should be pairs of files" << std::endl;
            return false;
        }
        else if ( min_kmers.size() > 0 && max_error.size() > 0 )
        {
            std::cerr << "--min-kmers and --max-error are mutually exclusive, please use just one" << std::endl;
            return false;
        }
        else if ( paired_mode > 1 )
        {
            std::cerr << "invalid --paired-mode" << std::endl;
            return false;
        }

        if ( reads_paired.size() > 1 && paired_mode == -1 )
            paired_mode = 1;

        output_unclassified    = !output_unclassified_file.empty() ? true : false;
        uint16_t extra_threads = output_unclassified ? 1 : 0;
        if ( threads <= min_threads )
            threads = min_threads + extra_threads;
        clas_threads = threads - 2 - extra_threads; //-1 reading, -1 printing clasified, -1 printing unclassified

        if ( n_batches < 1 )
            n_batches = 1;

        if ( n_reads < 1 )
            n_reads = 1;

        // default min_kmers
        if ( min_kmers.size() == 0 && max_error.size() == 0 )
            min_kmers.push_back( default_min_kmers );

        return parse_hierarchy();
    }

    bool parse_hierarchy()
    {

        if ( bloom_filter_files.size() != group_bin_files.size() )
        {
            std::cerr << "Filter and group-bin files do not match" << std::endl;
            return false;
        }

        if ( filter_hierarchy.size() == 1 && bloom_filter_files.size() > 1 )
        {
            for ( uint16_t b = 1; b < bloom_filter_files.size(); ++b )
            {
                filter_hierarchy.push_back( filter_hierarchy[0] );
            }
        }
        else if ( filter_hierarchy.size() != bloom_filter_files.size() )
        {
            std::cerr << "Hierarchy does not match with the number of provided files" << std::endl;
            return false;
        }

        if ( min_kmers.size() > 0 )
        {
            // If only one max error was given, repeat it for every filter
            if ( min_kmers.size() == 1 && bloom_filter_files.size() > 1 )
            {
                for ( uint16_t b = 1; b < bloom_filter_files.size(); ++b )
                {
                    min_kmers.push_back( min_kmers[0] );
                }
            }
            else if ( min_kmers.size() != bloom_filter_files.size() )
            {
                std::cerr << "Please provide a single or one-per-filter --min-kmers value[s]" << std::endl;
                return false;
            }
        }
        else
        {
            // If only one max error was given, repeat it for every filter
            if ( max_error.size() == 1 && bloom_filter_files.size() > 1 )
            {
                for ( uint16_t b = 1; b < bloom_filter_files.size(); ++b )
                {
                    max_error.push_back( max_error[0] );
                }
            }
            else if ( max_error.size() != bloom_filter_files.size() )
            {
                std::cerr << "Please provide a single or one-per-filter --max-error value[s]" << std::endl;
                return false;
            }
        }

        // If only one max error was given, repeat it for every filter
        std::vector< std::string > sorted_hierarchy = filter_hierarchy;
        std::sort( sorted_hierarchy.begin(), sorted_hierarchy.end() );
        uint16_t unique_hierarchy =
            std::unique( sorted_hierarchy.begin(), sorted_hierarchy.end() ) - sorted_hierarchy.begin();

        if ( max_error_unique.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                max_error_unique.push_back( max_error_unique[0] );
            }
        }
        else if ( max_error_unique.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --max-error-unique value[s]" << std::endl;
            return false;
        }

        uint16_t hierarchy_count = 0;
        for ( uint16_t h = 0; h < filter_hierarchy.size(); ++h )
        {
            auto filter_cfg = FilterConfig{ bloom_filter_files[h],
                                            group_bin_files[h],
                                            static_cast<int16_t>( max_error.size() > 0 ? max_error[h] : -1 ),
                                            ( min_kmers.size() > 0 ? min_kmers[h] : -1 ) };

            if ( h_filters.find( filter_hierarchy[h] ) == h_filters.end() )
            { // not found
                std::vector< FilterConfig > fc;
                fc.push_back( filter_cfg );

                std::string final_output_file;
                if ( !split_output_file_hierarchy || ( !output_file.empty() && unique_hierarchy == 1 ) )
                    final_output_file = output_file;
                else if ( !output_file.empty() && unique_hierarchy > 1 )
                    final_output_file = output_file + "_" + filter_hierarchy[h];
                else
                    final_output_file = "";

                h_filters[filter_hierarchy[h]] =
                    HierarchyConfig{ fc, max_error_unique[hierarchy_count], final_output_file };

                ++hierarchy_count;
            }
            else
            { // found
                h_filters[filter_hierarchy[h]].filters.push_back( filter_cfg );
            }
        }
        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };

    stream << newl;
    for ( auto const& hierarchy_config : config.h_filters )
    {
        stream << hierarchy_config.first << ") max-error-unique: " << hierarchy_config.second.max_error_unique
               << ", output-file: " << hierarchy_config.second.output_file << newl;
        for ( auto const& filter_config : hierarchy_config.second.filters )
        {
            stream << "  bloom-filter: " << filter_config.bloom_filter_file
                   << ", group-bin: " << filter_config.group_bin_file << ", min-kmers: " << filter_config.min_kmers
                   << ", max-error: " << filter_config.max_error << newl;
        }
    }
    stream << newl;
    stream << "--n-batches                 " << config.n_batches << newl;
    stream << "--n-reads                   " << config.n_reads << newl;
    stream << "--offset                    " << config.offset << newl;
    stream << "--verbose                   " << config.verbose << newl;
    stream << "--threads                   " << config.threads << newl;
    stream << "--paired-mode               " << config.paired_mode << newl;
    stream << "--output-unclassified-file  " << config.output_unclassified_file << newl;
    stream << "--reads-single              " << newl;
    for ( const auto& s : config.reads_single )
        stream << "                            " << s << newl;
    stream << "--reads-paired                  " << newl;
    for ( const auto& s : config.reads_paired )
        stream << "                            " << s << newl;
    stream << newl;

    return stream;
}

} // namespace GanonClassify
