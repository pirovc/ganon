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
    FilterConfig()
    {
    }

    FilterConfig( std::string _ibf_file, int16_t _max_error, float _min_kmers )
    {
        ibf_file  = _ibf_file;
        max_error = _max_error;
        min_kmers = _min_kmers;
    }

    std::string ibf_file;
    std::string map_file = "";
    std::string tax_file = "";
    int16_t     max_error;
    float       min_kmers;
};

struct HierarchyConfig
{
    std::vector< FilterConfig > filters;
    int16_t                     max_error_unique;
    uint8_t                     kmer_size;
    int16_t                     strata_filter;
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

    std::vector< std::string > hierarchy_labels{ "1_default" };

    std::vector< uint8_t > kmer_size{ 19 };
    std::vector< float >   min_kmers{ 0.25 };
    std::vector< int16_t > max_error;
    std::vector< int16_t > max_error_unique{ -1 };
    std::vector< int16_t > strata_filter{ 0 };
    uint8_t                offset = 1;

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

    bool validate()
    {
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
        for ( uint16_t i = 0; i < min_kmers.size(); ++i )
        {
            if ( min_kmers[i] < 0 || min_kmers[i] > 1 )
            {
                valid_val = false;
                break;
            }
        }
        if ( !valid_val )
        {
            std::cerr << "--min-kmers value should be between 0 and 1" << std::endl;
            return false;
        }

        valid_val = true;
        for ( uint16_t i = 0; i < max_error.size(); ++i )
        {
            if ( max_error[i] < 0 )
            {
                valid_val = false;
                break;
            }
        }
        if ( !valid_val )
        {
            std::cerr << "--max-error value should be greater than 0" << std::endl;
            return false;
        }

        // default is min_kmers, if min_error is provided, use it
        if ( max_error.size() > 0 )
            min_kmers[0] = -1; // reset min_kmers
        else
            max_error.push_back( -1 ); // reset max_error

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

        if ( offset < 1 )
            offset = 1;

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
        if ( max_error.size() == 1 && ibf.size() > 1 )
        {
            for ( uint16_t b = 1; b < ibf.size(); ++b )
            {
                max_error.push_back( max_error[0] );
            }
        }
        else if ( max_error.size() != ibf.size() )
        {
            std::cerr << "Please provide a single or one-per-filter --max-error value[s]" << std::endl;
            return false;
        }

        // If only one max error was given, repeat it for every filter
        if ( min_kmers.size() == 1 && ibf.size() > 1 )
        {
            for ( uint16_t b = 1; b < ibf.size(); ++b )
            {
                min_kmers.push_back( min_kmers[0] );
            }
        }
        else if ( min_kmers.size() != ibf.size() )
        {
            std::cerr << "Please provide a single or one-per-filter --min-kmers value[s]" << std::endl;
            return false;
        }


        std::vector< std::string > sorted_hierarchy = hierarchy_labels;
        std::sort( sorted_hierarchy.begin(), sorted_hierarchy.end() );
        // get unique hierarcy labels
        uint16_t unique_hierarchy =
            std::unique( sorted_hierarchy.begin(), sorted_hierarchy.end() ) - sorted_hierarchy.begin();

        for ( uint16_t k = 0; k < kmer_size.size(); ++k )
        {
            if ( offset > kmer_size[k] )
            {
                std::cerr << "Offset cannot be bigger than k-mer size" << std::endl;
                return false;
            }
        }
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

        if ( strata_filter.size() == 1 && unique_hierarchy > 1 )
        {
            for ( uint16_t b = 1; b < unique_hierarchy; ++b )
            {
                strata_filter.push_back( strata_filter[0] );
            }
        }
        else if ( strata_filter.size() != unique_hierarchy )
        {
            std::cerr << "Please provide a single or one-per-hierarchy --strata-filter value[s]" << std::endl;
            return false;
        }

        uint16_t hierarchy_count = 0;
        for ( uint16_t h = 0; h < hierarchy_labels.size(); ++h )
        {
            auto filter_cfg = FilterConfig{ ibf[h], max_error[h], min_kmers[h] };
            if ( map.size() > 0 )
                filter_cfg.map_file = map[h];
            if ( tax.size() > 0 )
                filter_cfg.tax_file = tax[h];

            if ( parsed_hierarchy.find( hierarchy_labels[h] ) == parsed_hierarchy.end() )
            { // not found
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
                                                                         max_error_unique[hierarchy_count],
                                                                         kmer_size[hierarchy_count],
                                                                         strata_filter[hierarchy_count],
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
            stream << " --kmer-size: " << hierarchy_config.second.kmer_size << newl;
            stream << " --max-error-unique: " << hierarchy_config.second.max_error_unique << newl;
            stream << " --strata-filter: " << hierarchy_config.second.strata_filter << newl;
        }
        for ( auto const& filter_config : hierarchy_config.second.filters )
        {
            if ( filter_config.min_kmers > -1 )
                stream << "  --min-kmers: " << filter_config.min_kmers << newl;
            if ( filter_config.max_error > -1 )
                stream << "  --max-error: " << filter_config.max_error << newl;
            stream << "  --ibf: " << filter_config.ibf_file << newl;
            stream << "  --map: " << filter_config.map_file << newl;
            stream << "  --tax: " << filter_config.tax_file << newl;
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
    stream << "--reads-single              " << newl;
    for ( const auto& s : config.single_reads )
        stream << "                            " << s << newl;
    stream << "--reads-paired                  " << newl;
    for ( const auto& s : config.paired_reads )
        stream << "                            " << s << newl;
    stream << newl;
    stream << "Parameters:" << newl;
    stream << "--offset                    " << unsigned( config.offset ) << newl;
    stream << "--output-prefix             " << config.output_prefix << newl;
    stream << "--output-lca                " << config.output_lca << newl;
    stream << "--output-all                " << config.output_all << newl;
    stream << "--output-unclassified       " << config.output_unclassified << newl;
    stream << "--output-single             " << config.output_single << newl;
    stream << "--threads                   " << config.threads << newl;
    stream << "--n-batches                 " << config.n_batches << newl;
    stream << "--n-reads                   " << config.n_reads << newl;
    stream << "--verbose                   " << config.verbose << newl;
    stream << "--quiet                     " << config.quiet << newl;
    stream << separator << newl;

    return stream;
}

} // namespace GanonClassify
