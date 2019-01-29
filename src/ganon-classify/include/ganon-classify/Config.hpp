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
    uint16_t    max_error;
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
    std::vector< std::string > bloom_filter_files;
    std::vector< std::string > group_bin_files;
    std::vector< std::string > reads;

    std::string max_error                   = "3";
    std::string max_error_unique            = "-1";
    std::string output_file                 = "";
    std::string output_unclassified_file    = "";
    std::string filter_hierarchy            = "1";
    uint16_t    offset                      = 1;
    uint16_t    threads                     = 3;
    bool        verbose                     = false;
    bool        split_output_file_hierarchy = false;

    // private:
    const uint16_t                           min_threads = 3;
    uint16_t                                 clas_threads;
    bool                                     output_unclassified = false;
    bool                                     unique_filtering    = false;
    std::map< std::string, HierarchyConfig > h_filters;


    bool validate()
    {
        output_unclassified    = !output_unclassified_file.empty() ? true : false;
        uint16_t extra_threads = output_unclassified ? 1 : 0;
        if ( threads <= min_threads )
            threads = min_threads + extra_threads;
        clas_threads = threads - 2 - extra_threads; //-1 reading, -1 printing clasified, -1 printing unclassified

        return parse_hierarchy();
    }

    std::vector< std::string > split( const std::string& s, char delimiter )
    {
        std::vector< std::string > tokens;
        std::string                token;
        std::istringstream         tokenStream( s );
        while ( std::getline( tokenStream, token, delimiter ) )
            tokens.push_back( token );
        return tokens;
    }


    bool parse_hierarchy()
    {

        if ( bloom_filter_files.size() != group_bin_files.size() )
        {
            std::cerr << "Filter and group-bin files do not match" << std::endl;
            return false;
        }

        std::vector< std::string > max_errors        = split( max_error, ',' );
        std::vector< std::string > max_errors_unique = split( max_error_unique, ',' );
        std::vector< std::string > hierarchy         = split( filter_hierarchy, ',' );

        if ( hierarchy.size() == 1 && bloom_filter_files.size() > 1 )
        {
            for ( uint16_t b = 1; b < bloom_filter_files.size(); ++b )
            {
                hierarchy.push_back( hierarchy[0] );
            }
        }
        else if ( hierarchy.size() != bloom_filter_files.size() )
        {
            std::cerr << "Hierarchy does not match with the number of provided files" << std::endl;
            return false;
        }

        // If only one max error was given, repeat it for every filter
        if ( max_errors.size() == 1 && bloom_filter_files.size() > 1 )
        {
            for ( uint16_t b = 1; b < bloom_filter_files.size(); ++b )
            {
                max_errors.push_back( max_errors[0] );
            }
        }
        else if ( max_errors.size() != bloom_filter_files.size() )
        {
            std::cerr << "Please give a single max error value or one-per-filter value" << std::endl;
            return false;
        }

        // If only one max error was given, repeat it for every filter
        std::vector< std::string > sorted_hierarchy = hierarchy;
        std::sort( sorted_hierarchy.begin(), sorted_hierarchy.end() );
        uint16_t unique_hierarchy =
            std::unique( sorted_hierarchy.begin(), sorted_hierarchy.end() ) - sorted_hierarchy.begin();

        if ( max_errors_unique.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                max_errors_unique.push_back( max_errors_unique[0] );
            }
        }
        else if ( max_errors_unique.size() != unique_hierarchy )
        {
            std::cerr << "Please give a single max error unique value or one-per-hierarchy level value" << std::endl;
            return false;
        }

        uint16_t hierarchy_count = 0;
        for ( uint16_t h = 0; h < hierarchy.size(); ++h )
        {

            if ( h_filters.find( hierarchy[h] ) == h_filters.end() )
            { // not found
                std::cerr << hierarchy[h] << "," << max_errors[h] << "," << max_errors_unique[hierarchy_count]
                          << std::endl;

                std::vector< FilterConfig > fc;
                fc.push_back( FilterConfig{ bloom_filter_files[h],
                                            group_bin_files[h],
                                            static_cast< uint16_t >( std::stoi( max_errors[h] ) ) } );

                std::string final_output_file;
                if ( !split_output_file_hierarchy || ( !output_file.empty() && unique_hierarchy == 1 ) )
                    final_output_file = output_file;
                else if ( !output_file.empty() && unique_hierarchy > 1 )
                    final_output_file = output_file + "_" + hierarchy[h];
                else
                    final_output_file = "";

                h_filters[hierarchy[h]] = HierarchyConfig{
                    fc, static_cast< int16_t >( std::stoi( max_errors_unique[hierarchy_count] ) ), final_output_file
                };

                ++hierarchy_count;
            }
            else
            { // found
                h_filters[hierarchy[h]].filters.push_back( FilterConfig{
                    bloom_filter_files[h], group_bin_files[h], static_cast< uint16_t >( std::stoi( max_errors[h] ) ) } );
            }
        }
        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };

    stream << newl;
    stream << "--bloom-filter          " << newl;
    for ( const auto& s : config.bloom_filter_files )
        stream << "                            " << s << newl;
    stream << "--group-bin          " << newl;
    for ( const auto& s : config.group_bin_files )
        stream << "                            " << s << newl;
    stream << "--filter-hierarchy          " << config.filter_hierarchy << newl;
    stream << "--max-error                 " << config.max_error << newl;
    stream << "--max-error-unique          " << config.max_error_unique << newl;
    stream << "--offset                    " << config.offset << newl;
    stream << "--output-file               " << config.output_file << newl;
    stream << "--output-unclassified-file  " << config.output_unclassified_file << newl;
    stream << "--verbose                   " << config.verbose << newl;
    stream << "--threads                   " << config.threads << newl;
    stream << "--reads                     " << newl;
    for ( const auto& s : config.reads )
        stream << "                            " << s << newl;
    stream << newl;
    stream << "(hierarchy)  " << newl;
    for ( auto const& hierarchy_config : config.h_filters )
    {
        stream << hierarchy_config.first << ") max-error-unique: " << hierarchy_config.second.max_error_unique
               << ", output-file: " << hierarchy_config.second.output_file << newl;
        for ( auto const& filter_config : hierarchy_config.second.filters )
        {
            stream << "  bloom-filter: " << filter_config.bloom_filter_file
                   << ", group-bin: " << filter_config.group_bin_file << ", max-error: " << filter_config.max_error
                   << newl;
        }
    }
    stream << newl;

    return stream;
}

} // namespace GanonClassify
