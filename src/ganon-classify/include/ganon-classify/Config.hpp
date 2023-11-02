#pragma once

#include <algorithm>
#include <cinttypes>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <seqan3/core/debug_stream.hpp>
#include <string>
#include <vector>

namespace GanonClassify
{


struct Config
{

public:
    std::vector< std::string > single_reads;
    std::vector< std::string > paired_reads;

    std::vector< std::string > ibf;
    std::vector< std::string > tax;

    std::vector< std::string > hierarchy_labels{ "H1" };

    std::vector< double > rel_cutoff{ 0.2 };
    std::vector< double > rel_filter{ 0.0 };
    std::vector< double > fpr_query{ 1.0 };

    std::string output_prefix       = "";
    bool        output_lca          = false;
    bool        output_all          = false;
    bool        output_unclassified = false;
    bool        output_single       = false;

    bool        hibf          = false;
    bool        skip_lca      = false;
    std::string tax_root_node = "1";
    uint16_t    threads       = 1;
    size_t      n_batches     = 1000;
    size_t      n_reads       = 400;
    bool        verbose       = false;
    bool        quiet         = false;


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

        valid_val = true;
        for ( uint16_t i = 0; i < fpr_query.size(); ++i )
        {
            if ( fpr_query[i] < 0 || fpr_query[i] > 1 )
            {
                valid_val = false;
                break;
            }
        }
        if ( !valid_val )
        {
            std::cerr << "--fpr-query values should be set between 0 and 1 (1 to disable)" << std::endl;
            return false;
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

        if ( validate_hierarchy() == false )
            return false;

        // Disable LCA without tax files
        if ( tax.size() == 0 )
            skip_lca = true;

        return true;
    }

    bool validate_hierarchy()
    {

        std::vector< std::string > sorted_hierarchy = hierarchy_labels;
        std::sort( sorted_hierarchy.begin(), sorted_hierarchy.end() );
        // get unique hierarcy labels
        uint16_t unique_hierarchy =
            std::unique( sorted_hierarchy.begin(), sorted_hierarchy.end() ) - sorted_hierarchy.begin();

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

        if ( fpr_query.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                fpr_query.push_back( fpr_query[0] );
            }
        }
        else if ( fpr_query.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --fpr-query value[s]" << std::endl;
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
            std::cerr << "--hierarchy does not match with the number of --ibf and --tax" << std::endl;
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

        return true;
    }
};

inline std::ostream& operator<<( std::ostream& stream, const Config& config )
{
    constexpr auto newl{ "\n" };
    constexpr auto separator{ "----------------------------------------------------------------------" };

    stream << separator << newl;
    if ( config.single_reads.size() )
    {
        stream << "--single-reads        " << newl;
        for ( const auto& s : config.single_reads )
            stream << "                      " << s << newl;
    }
    if ( config.paired_reads.size() )
    {
        stream << "--paired-reads        " << newl;
        for ( const auto& s : config.paired_reads )
            stream << "                      " << s << newl;
    }
    stream << "--output-prefix       " << config.output_prefix << newl;
    stream << "--output-lca          " << config.output_lca << newl;
    stream << "--output-all          " << config.output_all << newl;
    stream << "--output-unclassified " << config.output_unclassified << newl;
    stream << "--output-single       " << config.output_single << newl;
    stream << "--hibf                " << config.hibf << newl;
    stream << "--threads             " << config.threads << newl;
    stream << "--n-batches           " << config.n_batches << newl;
    stream << "--n-reads             " << config.n_reads << newl;
    stream << "--skip-lca            " << config.skip_lca << newl;
    stream << "--verbose             " << config.verbose << newl;
    stream << "--quiet               " << config.quiet << newl;
    stream << separator << newl;

    return stream;
}

} // namespace GanonClassify
