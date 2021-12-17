#pragma once

#include <algorithm>
#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>
#include <string>
#include <vector>

namespace GanonClassify
{

struct FilterConfig
{
    FilterConfig()
    {
    }

    FilterConfig( std::string _ibf_file, int16_t _abs_cutoff, double _rel_cutoff )
    {
        ibf_file   = _ibf_file;
        abs_cutoff = _abs_cutoff;
        rel_cutoff = _rel_cutoff;
    }

    std::string ibf_file;
    std::string map_file = "";
    std::string tax_file = "";
    int16_t     abs_cutoff;
    double      rel_cutoff;
};

struct HierarchyConfig
{
    std::vector< FilterConfig > filters;
    uint8_t                     kmer_size;
    uint8_t                     window_size;
    uint8_t                     offset;
    double                      rel_filter;
    int16_t                     abs_filter;
    std::string                 output_file_lca;
    std::string                 output_file_all;
};

struct Config
{

public:
    std::vector< std::string > single_reads;
    std::vector< std::string > paired_reads;

    std::vector< std::string > ibf;
    std::vector< std::string > map;
    std::vector< std::string > tax;

    std::vector< std::string > hierarchy_labels{ "H1" };

    std::vector< uint8_t > kmer_size{ 19 };
    std::vector< uint8_t > window_size{ 0 };
    std::vector< uint8_t > offset{ 1 };

    std::vector< double >  rel_cutoff{ 0.25 };
    std::vector< int16_t > abs_cutoff;
    std::vector< double >  rel_filter;
    std::vector< int16_t > abs_filter{ 0 };

    std::string output_prefix       = "";
    bool        output_lca          = false;
    bool        output_all          = false;
    bool        output_unclassified = false;
    bool        output_single       = false;

    uint16_t threads   = 3;
    uint32_t n_batches = 1000;
    uint32_t n_reads   = 400;
    bool     verbose   = false;
    bool     quiet     = false;

    uint16_t                                 threads_classify;
    std::map< std::string, HierarchyConfig > parsed_hierarchy;


    bool check_files( std::vector< std::string > const& files )
    {
        for ( auto const& file : files )
        {
            if ( !std::filesystem::exists( file ) )
            {
                std::cerr << "file not found: " << file << std::endl;
                return false;
            }
            else if ( std::filesystem::file_size( file ) == 0 )
            {
                std::cerr << "file is empty: " << file << std::endl;
                return false;
            }
        }
        return true;
    }


    bool validate()
    {
        if ( !check_files( single_reads ) )
            return false;
        if ( !check_files( paired_reads ) )
            return false;
        if ( !check_files( ibf ) )
            return false;
        if ( !check_files( map ) )
            return false;
        if ( !check_files( tax ) )
            return false;

        if ( ibf.size() == 0 || ( paired_reads.size() == 0 && single_reads.size() == 0 ) )
        {
            std::cerr << "--ibf and --[single|paired]-reads are mandatory" << std::endl;
            return false;
        }

        if ( paired_reads.size() % 2 != 0 )
        {
            std::cerr << "--paired-reads should be an even number of files (pairs)" << std::endl;
            return false;
        }

        if ( tax.size() > 0 && map.size() == 0 )
        {
            std::cerr << "--map is required to use --tax" << std::endl;
            return false;
        }

        bool valid_val = true;
        for ( uint16_t i = 0; i < rel_cutoff.size(); ++i )
        {
            if ( rel_cutoff[i] < 0 || rel_cutoff[i] > 1 )
            {
                valid_val = false;
                break;
            }
        }
        if ( !valid_val )
        {
            std::cerr << "--rel-cutoff values should be set between 0 and 1 (0 to disable)" << std::endl;
            return false;
        }

        valid_val = true;
        for ( uint16_t i = 0; i < abs_filter.size(); ++i )
        {
            if ( abs_filter[i] < 0 && abs_filter[i] != -1 )
            {
                valid_val = false;
                break;
            }
        }
        if ( !valid_val )
        {
            std::cerr << "--abs-filter values should be >= 0 (-1 to disable)" << std::endl;
            return false;
        }

        valid_val = true;
        for ( uint16_t i = 0; i < abs_cutoff.size(); ++i )
        {
            if ( abs_cutoff[i] < 0 && abs_cutoff[i] != -1 )
            {
                valid_val = false;
                break;
            }
        }
        if ( !valid_val )
        {
            std::cerr << "--abs-cutoff values should be >= 0 (-1 to disable)" << std::endl;
            return false;
        }

        valid_val = true;
        for ( uint16_t i = 0; i < rel_filter.size(); ++i )
        {
            if ( rel_filter[i] < 0 || rel_filter[i] > 1 )
            {
                valid_val = false;
                break;
            }
        }
        if ( !valid_val )
        {
            std::cerr << "--rel-filter values should be set between 0 and 1 (1 to disable)" << std::endl;
            return false;
        }

        // default is rel_cutoff
        if ( abs_cutoff.size() > 0 )
            rel_cutoff[0] = -1; // reset rel_cutoff
        else
            abs_cutoff.push_back( -1 ); // reset abs_cutoff

        // default is abs_filter
        if ( rel_filter.size() > 0 )
            abs_filter[0] = -1; // reset abs_filter
        else
            rel_filter.push_back( -1 ); // reset rel_filter

        if ( threads <= 3 )
        {
            threads_classify = 1;
        }
        else
        {
            threads_classify = threads - 2;
        }

        if ( n_batches < 1 )
            n_batches = 1;

        if ( n_reads < 1 )
            n_reads = 1;

        if ( output_prefix.empty() )
        {
            output_lca          = false;
            output_all          = false;
            output_unclassified = false;
        }

        return parse_hierarchy();
    }

    bool parse_hierarchy()
    {
        if ( map.size() > 0 && ibf.size() != map.size() )
        {
            std::cerr << "The number of files provided with --ibf and --map should match" << std::endl;
            return false;
        }

        if ( tax.size() > 0 && ibf.size() != tax.size() )
        {
            std::cerr << "The number of files provided with --ibf and --tax should match" << std::endl;
            return false;
        }

        // One hierarchy value for multiple hierarchies
        if ( hierarchy_labels.size() == 1 && ibf.size() > 1 )
        {
            for ( uint16_t b = 1; b < ibf.size(); ++b )
            {
                hierarchy_labels.push_back( hierarchy_labels[0] );
            }
        }
        else if ( hierarchy_labels.size() != ibf.size() )
        {
            std::cerr << "--hierarchy does not match with the number of --ibf, --map and --tax" << std::endl;
            return false;
        }

        // If only one max error was given, repeat it for every filter
        if ( abs_cutoff.size() == 1 && ibf.size() > 1 )
        {
            for ( uint16_t b = 1; b < ibf.size(); ++b )
            {
                abs_cutoff.push_back( abs_cutoff[0] );
            }
        }
        else if ( abs_cutoff.size() != ibf.size() )
        {
            std::cerr << "Please provide a single or one-per-filter --abs-cutoff value[s]" << std::endl;
            return false;
        }

        // If only one max error was given, repeat it for every filter
        if ( rel_cutoff.size() == 1 && ibf.size() > 1 )
        {
            for ( uint16_t b = 1; b < ibf.size(); ++b )
            {
                rel_cutoff.push_back( rel_cutoff[0] );
            }
        }
        else if ( rel_cutoff.size() != ibf.size() )
        {
            std::cerr << "Please provide a single or one-per-filter --rel-cutoff value[s]" << std::endl;
            return false;
        }

        std::vector< std::string > sorted_hierarchy = hierarchy_labels;
        std::sort( sorted_hierarchy.begin(), sorted_hierarchy.end() );
        // get unique hierarcy labels
        uint16_t unique_hierarchy =
            std::unique( sorted_hierarchy.begin(), sorted_hierarchy.end() ) - sorted_hierarchy.begin();

        if ( kmer_size.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                kmer_size.push_back( kmer_size[0] );
            }
        }
        else if ( kmer_size.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --kmer-size value[s]" << std::endl;
            return false;
        }

        for ( uint16_t w = 0; w < window_size.size(); ++w )
        {
            if ( window_size[w] > 0 && window_size[w] < kmer_size[w] )
            {
                std::cerr << "--window-size has to be >= --kmer-size (or 0 to disable windows)" << std::endl;
                return false;
            }
        }
        if ( window_size.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                window_size.push_back( window_size[0] );
            }
        }
        else if ( window_size.size() > 0 && window_size.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --window-size value[s]" << std::endl;
            return false;
        }

        if ( offset.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                offset.push_back( offset[0] );
            }
        }
        else if ( offset.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --offset value[s]" << std::endl;
            return false;
        }

        if ( rel_filter.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                rel_filter.push_back( rel_filter[0] );
            }
        }
        else if ( rel_filter.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --rel-filter value[s]" << std::endl;
            return false;
        }

        if ( abs_filter.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                abs_filter.push_back( abs_filter[0] );
            }
        }
        else if ( abs_filter.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --abs-filter value[s]" << std::endl;
            return false;
        }

        uint16_t hierarchy_count = 0;
        for ( uint16_t h = 0; h < hierarchy_labels.size(); ++h )
        {

            auto filter_cfg = FilterConfig{ ibf[h], abs_cutoff[h], rel_cutoff[h] };
            if ( map.size() > 0 )
                filter_cfg.map_file = map[h];
            if ( tax.size() > 0 )
                filter_cfg.tax_file = tax[h];

            if ( parsed_hierarchy.find( hierarchy_labels[h] ) == parsed_hierarchy.end() )
            { // not found
                // validate by hiearchy
                // if window size is used in this level
                if ( window_size[hierarchy_count] > 0 )
                {
                    if ( abs_filter[hierarchy_count] != -1 || abs_cutoff[h] != -1 )
                    {
                        std::cerr << "minimizers (--window-size) can only be used with --rel-cutoff and --rel-filter"
                                  << std::endl;
                        return false;
                    }

                    if ( offset[hierarchy_count] > 1 )
                    {
                        std::cerr << "minimizers (--window-size) can not be used with --offset" << std::endl;
                        return false;
                    }
                }

                std::vector< FilterConfig > fc;
                fc.push_back( filter_cfg );
                std::string output_file_lca = "";
                std::string output_file_all = "";
                if ( !output_prefix.empty() && unique_hierarchy > 1 && !output_single )
                {
                    output_file_lca = output_prefix + "." + hierarchy_labels[h] + ".lca";
                    output_file_all = output_prefix + "." + hierarchy_labels[h] + ".all";
                }
                else if ( !output_prefix.empty() )
                {
                    output_file_lca = output_prefix + ".lca";
                    output_file_all = output_prefix + ".all";
                }

                parsed_hierarchy[hierarchy_labels[h]] = HierarchyConfig{ fc,
                                                                         kmer_size[hierarchy_count],
                                                                         window_size[hierarchy_count],
                                                                         offset[hierarchy_count],
                                                                         rel_filter[hierarchy_count],
                                                                         abs_filter[hierarchy_count],
                                                                         output_file_lca,
                                                                         output_file_all };
                ++hierarchy_count;
            }
            else
            { // found
                parsed_hierarchy[hierarchy_labels[h]].filters.push_back( filter_cfg );
            }
        }
        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };
    constexpr auto separator{ "----------------------------------------------------------------------" };

    stream << separator << newl;
    stream << "Database hierarchy:" << newl;
    for ( auto const& hierarchy_config : config.parsed_hierarchy )
    {
        if ( !hierarchy_config.first.empty() )
        {
            stream << hierarchy_config.first << ")" << newl;
            stream << " --kmer-size:         " << unsigned( hierarchy_config.second.kmer_size ) << newl;
            if ( hierarchy_config.second.offset > 1 )
                stream << " --offset:            " << unsigned( hierarchy_config.second.offset ) << newl;
            if ( hierarchy_config.second.window_size > 0 )
                stream << " --window-size:       " << unsigned( hierarchy_config.second.window_size ) << newl;
            if ( hierarchy_config.second.rel_filter > -1 )
                stream << " --rel-filter:        " << hierarchy_config.second.rel_filter << newl;
            if ( hierarchy_config.second.abs_filter > -1 )
                stream << " --abs-filter:        " << hierarchy_config.second.abs_filter << newl;
        }
        for ( auto const& filter_config : hierarchy_config.second.filters )
        {
            stream << "   --ibf:             " << filter_config.ibf_file << newl;
            if ( !filter_config.map_file.empty() )
                stream << "   --map:             " << filter_config.map_file << newl;
            if ( !filter_config.tax_file.empty() )
                stream << "   --tax:             " << filter_config.tax_file << newl;
            if ( filter_config.rel_cutoff > -1 )
                stream << "   --rel-cutoff:      " << filter_config.rel_cutoff << newl;
            if ( filter_config.abs_cutoff > -1 )
                stream << "   --abs-cutoff:      " << filter_config.abs_cutoff << newl;
        }
        if ( !config.output_prefix.empty() )
        {
            stream << "  Output files:" << newl;
            stream << "    " << config.output_prefix + ".rep" << newl;
            if ( config.output_lca )
                stream << "    " << hierarchy_config.second.output_file_lca << newl;
            if ( config.output_all )
                stream << "    " << hierarchy_config.second.output_file_all << newl;
        }
    }
    stream << newl;
    stream << "Input files:" << newl;
    if ( config.single_reads.size() )
    {
        stream << "--reads-single        " << newl;
        for ( const auto& s : config.single_reads )
            stream << "                      " << s << newl;
    }
    if ( config.paired_reads.size() )
    {
        stream << "--reads-paired        " << newl;
        for ( const auto& s : config.paired_reads )
            stream << "                      " << s << newl;
    }
    stream << newl;
    stream << "Parameters:" << newl;
    if ( config.output_prefix.size() )
        stream << "--output-prefix       " << config.output_prefix << newl;
    stream << "--output-lca          " << config.output_lca << newl;
    stream << "--output-all          " << config.output_all << newl;
    stream << "--output-unclassified " << config.output_unclassified << newl;
    stream << "--output-single       " << config.output_single << newl;
    stream << "--threads             " << config.threads << newl;
    stream << "--n-batches           " << config.n_batches << newl;
    stream << "--n-reads             " << config.n_reads << newl;
    stream << "--verbose             " << config.verbose << newl;
    stream << "--quiet               " << config.quiet << newl;
    stream << separator << newl;

    return stream;
}

} // namespace GanonClassify
