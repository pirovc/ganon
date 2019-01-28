#pragma once

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
    uint16_t                    max_error_unique;
    std::string                 output_file;
};

struct Config
{

public:
    std::vector< std::string > bloom_filter_files;
    std::vector< std::string > group_bin_files;
    std::string                max_error        = "3";
    std::string                max_error_unique = "-1";

    std::string output_file;
    std::string output_unclassified_file;
    std::string filter_hierarchy;

    uint16_t                   offset  = 1;
    uint16_t                   threads = 3;
    std::vector< std::string > reads;
    bool                       verbose = false;

    // private:
    const uint16_t                           min_threads = 3;
    uint16_t                                 clas_threads;
    bool                                     output_unclassified = false;
    bool                                     unique_filtering    = false;
    std::map< std::string, HierarchyConfig > h_filters;


    bool validate()
    {
        threads             = min_threads < min_threads ? min_threads : threads;
        threads             = !output_unclassified_file.empty() ? threads + 1u : threads;
        output_unclassified = !output_unclassified_file.empty() ? true : false;
        //-1 reading, -1 printing clasified, -1 printing unclassified
        clas_threads = !output_unclassified_file.empty() ? threads - 3u : threads - 2u;

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
        std::vector< std::string > max_errors        = split( max_error, ',' );
        std::vector< std::string > max_errors_unique = split( max_error_unique, ',' );

        if ( filter_hierarchy.empty() )
        {

            if ( bloom_filter_files.size() != group_bin_files.size() )
            {
                std::cerr << "Filter and group-bin files do not match" << std::endl;
                return false;
            }
            else
            {
                std::vector< FilterConfig > fc;
                for ( uint16_t h = 0; h < bloom_filter_files.size(); ++h )
                {
                    fc.push_back(
                        FilterConfig{ bloom_filter_files[h], group_bin_files[h], std::stoi( max_errors[h] ) } );
                }
                h_filters["1"] = HierarchyConfig{ fc, std::stoi( max_errors_unique[0] ), output_file };
            }
        }
        else
        {
            std::vector< std::string > hierarchy = split( filter_hierarchy, ',' );

            if ( hierarchy.size() != bloom_filter_files.size() || hierarchy.size() != group_bin_files.size() )
            {
                std::cerr << "Hierarchy does not match with the number of provided files" << std::endl;
                return false;
            }
            else
            {
                for ( uint16_t h = 0; h < hierarchy.size(); ++h )

                    if ( h_filters.find( hierarchy[h] ) == h_filters.end() )
                    { // not found
                        std::vector< FilterConfig > fc;
                        fc.push_back(
                            FilterConfig{ bloom_filter_files[h], group_bin_files[h], std::stoi( max_errors[h] ) } );

                        h_filters[hierarchy[h]] =
                            HierarchyConfig{ fc, std::stoi( max_errors_unique[0] ), hierarchy[h] + "_" + output_file };
                    }
                    else
                    { // found
                        h_filters[hierarchy[h]].filters.push_back(
                            FilterConfig{ bloom_filter_files[h], group_bin_files[h], std::stoi( max_errors[h] ) } );
                    }
            }
        }
        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };

    stream << newl;
    stream << "--filter-hierarchy          " << config.filter_hierarchy << newl;
    stream << "--bloom-filter,--group-bin  " << newl;
    for ( auto const& hierarchy_config : config.h_filters )
    {
        for ( auto const& filter_config : hierarchy_config.second.filters )
        {
            stream << "                            (" << hierarchy_config.first << ") "
                   << filter_config.bloom_filter_file << ", " << filter_config.group_bin_file << newl;
        }
    }
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

    return stream;
}

} // namespace GanonClassify
