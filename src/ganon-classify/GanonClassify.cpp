#include "GanonClassify.hpp"
#include "hierarchical_interleaved_bloom_filter.hpp"

#include <robin_hood.h>

#include <utils/IBFConfig.hpp>
#include <utils/LCA.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>
#include <utils/adjust_seed.hpp>
#include <utils/dna4_traits.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <cinttypes>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>


namespace GanonClassify
{

namespace detail
{

#ifdef LONGREADS
typedef uint32_t TIntCount;
#else
typedef uint16_t TIntCount;
#endif

typedef raptor::hierarchical_interleaved_bloom_filter< seqan3::data_layout::uncompressed > THIBF;
typedef seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed >              TIBF;
typedef robin_hood::unordered_map< std::string, std::tuple< size_t, double > >             TMatches;
typedef std::vector< std::tuple< size_t, std::string > >                                   TBinMap;
typedef robin_hood::unordered_map< std::string, std::vector< size_t > >                    TMap;
typedef robin_hood::unordered_map< std::string, double >                                   TTargetFpr;
typedef std::map< std::string, std::vector< std::pair< std::string, std::string > > >      TReadConfig;

struct PairHash
{
    size_t operator()( const std::pair< std::string, std::string >& p ) const
    {
        size_t h1 = std::hash< std::string >{}( p.first );
        size_t h2 = std::hash< std::string >{}( p.second );
        return h1 ^ ( h2 << 1 );
    }
};

struct Node
{
    std::string parent;
    std::string rank;
    std::string name;
};


struct ReadBatches
{

    ReadBatches()
    {
    }

    ReadBatches( bool _paired )
    {
        paired = _paired;
    }

    ReadBatches( bool _paired, std::string _prefix )
    {
        paired = _paired;
        prefix = _prefix;
    }

    ReadBatches( bool _paired, std::vector< std::string > _ids, std::vector< std::vector< seqan3::dna4 > > _seqs )
    {
        paired = _paired;
        ids    = _ids;
        seqs   = _seqs;
    }

    ReadBatches( bool                                       _paired,
                 std::vector< std::string >                 _ids,
                 std::vector< std::vector< seqan3::dna4 > > _seqs,
                 std::vector< std::vector< seqan3::dna4 > > _seqs2 )
    {
        paired = _paired;
        ids    = _ids;
        seqs   = _seqs;
        seqs2  = _seqs2;
    }

    bool                                       paired = false;
    std::vector< std::string >                 ids;
    std::vector< std::vector< seqan3::dna4 > > seqs;
    std::vector< std::vector< seqan3::dna4 > > seqs2{};
    std::string                                prefix;
};

struct ReadMatch
{
    ReadMatch()
    {
    }

    ReadMatch( std::string _target, size_t _kmer_count )
    {
        target     = _target;
        kmer_count = _kmer_count;
    }

    std::string target;
    size_t      kmer_count;
};

struct ReadOut
{
    ReadOut()
    {
    }

    ReadOut( std::string _readID, std::string _prefix )
    {
        readID = _readID;
        prefix = _prefix;
    }

    std::string              readID;
    std::vector< ReadMatch > matches;
    std::string              prefix;
};

struct Rep
{
    size_t matches                    = 0;
    size_t seqs_lca                   = 0;
    size_t seqs_unique                = 0;
    size_t discarded_matches_filter   = 0;
    size_t discarded_matches_fprquery = 0;
};

struct Total
{
    size_t input_seqs                 = 0;
    size_t seqs_processed             = 0;
    size_t seqs_skipped_big           = 0;
    size_t seqs_skipped_small         = 0;
    size_t length_processed           = 0;
    size_t kmers_processed            = 0;
    size_t seqs_classified            = 0;
    size_t kmers_matches              = 0;
    size_t kmers_from_classified_seqs = 0;
    size_t matches                    = 0;
    size_t seqs_unique                = 0;
    size_t discarded_matches_filter   = 0;
    size_t discarded_matches_fprquery = 0;
};


typedef robin_hood::unordered_map< std::pair< std::string, std::string >, Rep, PairHash > TRep;
typedef robin_hood::unordered_map< std::string, Total >                                   TTotal;
typedef robin_hood::unordered_map< std::string, Node >                                    TTax;


struct Stats
{
    // General total
    TTotal total;
    // Total for each hierarchy
    std::map< std::string, TTotal > hierarchy_total;

    // Total among many files
    size_t total_seqs_processed   = 0;
    size_t total_length_processed = 0;
    size_t total_kmers_processed  = 0;

    void add_totals( std::string hierarchy_label, std::vector< detail::TTotal > const& totals )
    {
        // add several totals (from threads) into the stats
        for ( auto const& total_thread : totals )
        {
            for ( auto const& [prefix, t] : total_thread )
            {
                total_seqs_processed += t.seqs_processed;
                total_length_processed += t.length_processed;
                total_kmers_processed += t.kmers_processed;

                total[prefix].seqs_processed += t.seqs_processed;
                total[prefix].seqs_skipped_big += t.seqs_skipped_big;
                total[prefix].seqs_skipped_small += t.seqs_skipped_small;
                total[prefix].length_processed += t.length_processed;
                total[prefix].kmers_processed += t.kmers_processed;
                total[prefix].seqs_classified += t.seqs_classified;
                total[prefix].kmers_matches += t.kmers_matches;
                total[prefix].kmers_from_classified_seqs += t.kmers_from_classified_seqs;

                hierarchy_total[hierarchy_label][prefix].seqs_processed += t.seqs_processed;
                hierarchy_total[hierarchy_label][prefix].seqs_skipped_big += t.seqs_skipped_big;
                hierarchy_total[hierarchy_label][prefix].seqs_skipped_small += t.seqs_skipped_small;
                hierarchy_total[hierarchy_label][prefix].length_processed += t.length_processed;
                hierarchy_total[hierarchy_label][prefix].kmers_processed += t.kmers_processed;
                hierarchy_total[hierarchy_label][prefix].seqs_classified += t.seqs_classified;
                hierarchy_total[hierarchy_label][prefix].kmers_matches += t.kmers_matches;
                hierarchy_total[hierarchy_label][prefix].kmers_from_classified_seqs += t.kmers_from_classified_seqs;
            }
        }
    }

    void add_reports( std::string hierarchy_label, TRep const& report )
    {
        // add values from reports to stats
        for ( auto const& [prefix_target, rep] : report )
        {
            auto prefix = prefix_target.first;

            total[prefix].matches += rep.matches;
            total[prefix].seqs_unique += rep.seqs_unique;
            total[prefix].discarded_matches_filter += rep.discarded_matches_filter;
            total[prefix].discarded_matches_fprquery += rep.discarded_matches_fprquery;

            hierarchy_total[hierarchy_label][prefix].matches += rep.matches;
            hierarchy_total[hierarchy_label][prefix].seqs_unique += rep.seqs_unique;
            hierarchy_total[hierarchy_label][prefix].discarded_matches_filter += rep.discarded_matches_filter;
            hierarchy_total[hierarchy_label][prefix].discarded_matches_fprquery += rep.discarded_matches_fprquery;
        }
    }
};

struct FilterConfig
{
    FilterConfig()
    {
    }

    FilterConfig( std::string _ibf_file, double _rel_cutoff )
    {
        ibf_file   = _ibf_file;
        rel_cutoff = _rel_cutoff;
    }

    std::string ibf_file;
    std::string tax_file = "";
    double      rel_cutoff;
    IBFConfig   ibf_config;
    TTargetFpr  target_fpr;
};

struct HierarchyConfig
{
    std::vector< FilterConfig > filters;
    uint8_t                     kmer_size;
    uint32_t                    window_size;
    double                      rel_filter;
    double                      fpr_query;
    std::string                 output_file_lca;
    std::string                 output_file_all;
};

template < typename TFilter >
struct Filter
{
    TFilter      ibf;
    TMap         map;
    TTax         tax;
    size_t       bin_count = 0;
    FilterConfig filter_config;
};

bool parse_reads_config( Config& config, TReadConfig& reads_config )
{


    if ( config.batch_reads.size() > 0 )
    {
        std::string line;
        for ( auto const& batch_file : config.batch_reads )
        {
            std::ifstream infile( batch_file );
            while ( std::getline( infile, line, '\n' ) )
            {
                std::istringstream         stream_line( line );
                std::vector< std::string > fields;
                std::string                field;
                while ( std::getline( stream_line, field, '\t' ) )
                {
                    fields.push_back( field ); // prefix <tab> file1 [<tab> file2]
                }

                if ( fields.size() <= 1 )
                {
                    std::cerr << "ERROR: invalid --batch-reads file (prefix <tab> file1 [<tab> file2])" << std::endl;
                    return false;
                }

                if ( fields.size() >= 2 )
                {
                    if ( !std::filesystem::exists( fields[1] ) || std::filesystem::file_size( fields[1] ) == 0 )
                    {
                        std::cerr << "ERROR: file not found/empty: " << fields[1] << std::endl;
                        return false;
                    }
                    if ( fields.size() == 3 )
                    {
                        if ( !std::filesystem::exists( fields[2] ) || std::filesystem::file_size( fields[2] ) == 0 )
                        {
                            std::cerr << "ERROR: file not found/empty: " << fields[2] << std::endl;
                            return false;
                        }
                        reads_config[fields[0]].push_back( { fields[1], fields[2] } );
                    }
                    else
                    {
                        reads_config[fields[0]].push_back( { fields[1], "" } );
                    }
                }
            }
        }
    }
    else
    {
        for ( auto const& reads_file : config.single_reads )
        {
            reads_config[""].push_back( { reads_file, "" } );
        }
        for ( size_t pair_cnt = 0; pair_cnt < config.paired_reads.size(); pair_cnt += 2 )
        {
            reads_config[""].push_back( { config.paired_reads[pair_cnt], config.paired_reads[pair_cnt + 1] } );
        }
    }
    return true;
}

std::map< std::string, HierarchyConfig > parse_hierarchy( Config& config )
{

    std::map< std::string, HierarchyConfig > parsed_hierarchy;

    std::vector< std::string > sorted_hierarchy = config.hierarchy_labels;
    std::sort( sorted_hierarchy.begin(), sorted_hierarchy.end() );
    // get unique hierarcy labels
    const size_t unique_hierarchy =
        std::unique( sorted_hierarchy.begin(), sorted_hierarchy.end() ) - sorted_hierarchy.begin();

    size_t hierarchy_count = 0;
    for ( size_t h = 0; h < config.hierarchy_labels.size(); ++h )
    {

        auto filter_cfg = FilterConfig{ config.ibf[h], config.rel_cutoff[h] };
        if ( config.tax.size() > 0 )
            filter_cfg.tax_file = config.tax[h];

        if ( parsed_hierarchy.find( config.hierarchy_labels[h] ) == parsed_hierarchy.end() )
        { // not found
            // validate by hiearchy
            std::vector< FilterConfig > fc;
            fc.push_back( filter_cfg );
            std::string output_file_lca = "one";
            std::string output_file_all = "all";
            if ( unique_hierarchy > 1 && !config.output_single )
            {
                output_file_lca = config.hierarchy_labels[h] + "." + output_file_lca;
                output_file_all = config.hierarchy_labels[h] + "." + output_file_all;
            }

            parsed_hierarchy[config.hierarchy_labels[h]] = HierarchyConfig{ fc,
                                                                            0,
                                                                            0,
                                                                            config.rel_filter[hierarchy_count],
                                                                            config.fpr_query[hierarchy_count],
                                                                            output_file_lca,
                                                                            output_file_all };
            ++hierarchy_count;
        }
        else
        { // found
            parsed_hierarchy[config.hierarchy_labels[h]].filters.push_back( filter_cfg );
        }
    }

    return parsed_hierarchy;
}

void print_hierarchy( auto const& parsed_hierarchy )
{
    constexpr auto newl{ "\n" };
    std::cerr << "Database(s):" << newl;
    for ( auto const& hierarchy_config : parsed_hierarchy )
    {
        std::cerr << hierarchy_config.first << ":" << newl;
        std::cerr << "--rel-filter " << hierarchy_config.second.rel_filter << newl;
        std::cerr << "--fpr-query " << hierarchy_config.second.fpr_query << newl;
        for ( auto const& filter_config : hierarchy_config.second.filters )
        {
            if ( filter_config.rel_cutoff > -1 )
                std::cerr << "--rel-cutoff " << filter_config.rel_cutoff;
            std::cerr << " " << filter_config.ibf_file;
            if ( !filter_config.tax_file.empty() )
                std::cerr << ", " << filter_config.tax_file;

            std::cerr << newl;
        }
    }
    std::cerr << "----------------------------------------------------------------------" << newl;
}

void print_reads_config( detail::TReadConfig const& reads_config )
{
    constexpr auto newl{ "\n" };
    std::cerr << "Sequence(s):" << newl;
    for ( auto const& [prefix, files] : reads_config )
    {
        if ( !prefix.empty() )
        {
            std::cerr << prefix << ":" << newl;
        }
        for ( auto const& [file1, file2] : files )
        {
            std::cerr << file1;
            if ( !file2.empty() )
                std::cerr << ", " << file2;
            std::cerr << newl;
        }
    }
    std::cerr << "----------------------------------------------------------------------" << newl;
}

void print_output_files( const Config&                                   config,
                         const std::map< std::string, HierarchyConfig >& parsed_hierarchy,
                         const TReadConfig&                              reads_config )
{
    constexpr auto newl{ "\n" };
    std::cerr << "Output file(s):" << newl;
    for ( auto& [prefix, files] : reads_config )
    {
        if ( !prefix.empty() )
        {
            std::cerr << prefix << ":" << newl;
        }
        std::cerr << config.output_prefix + prefix + ".rep" << newl;
        if ( config.output_unclassified )
            std::cerr << config.output_prefix + prefix + ".unc" << newl;
        if ( config.output_stats )
            std::cerr << config.output_prefix + prefix + ".sta" << newl;
        for ( auto& [hierarchy_label, hierarchy_config] : parsed_hierarchy )
        {
            if ( config.output_lca )
                std::cerr << config.output_prefix + prefix + "." + hierarchy_config.output_file_lca << newl;
            if ( config.output_all )
                std::cerr << config.output_prefix + prefix + "." + hierarchy_config.output_file_all << newl;
        }
    }
    std::cerr << "----------------------------------------------------------------------" << newl;
}

inline TRep sum_reports( std::vector< TRep > const& reports )
{
    detail::TRep report_sum;
    for ( auto const& report : reports )
    {
        for ( auto const& [prefix_target, rep] : report )
        {
            report_sum[prefix_target].matches += rep.matches;
            report_sum[prefix_target].seqs_lca += rep.seqs_lca;
            report_sum[prefix_target].seqs_unique += rep.seqs_unique;
            report_sum[prefix_target].discarded_matches_filter += rep.discarded_matches_filter;
            report_sum[prefix_target].discarded_matches_fprquery += rep.discarded_matches_fprquery;
        }
    }
    return report_sum;
}

inline size_t threshold_rel( size_t n_hashes, double p )
{
    return std::ceil( n_hashes * p );
}

// https://stackoverflow.com/questions/44718971/calculate-binomial-coffeficient-very-reliably
inline double binom( double n, double k ) noexcept
{
    return std::exp( std::lgamma( n + 1 ) - std::lgamma( n - k + 1 ) - std::lgamma( k + 1 ) );
}


void select_matches( Filter< TIBF >&        filter,
                     TMatches&              matches,
                     std::vector< size_t >& hashes,
                     auto&                  agent,
                     size_t                 threshold_cutoff,
                     size_t&                max_count_read,
                     size_t&                min_count_read,
                     size_t                 n_hashes )
{
    // Count every occurrence on IBF
    seqan3::counting_vector< detail::TIntCount > counts = agent.bulk_count( hashes );

    for ( auto const& [target, bins] : filter.map )
    {
        // Sum counts among bins (split target (user bins) into several technical bins)
        size_t summed_count = 0;
        for ( auto const& binno : bins )
        {
            summed_count += counts[binno];
        }
        // summed_count can be higher than n_hashes for matches in several technical bins
        if ( summed_count > n_hashes )
            summed_count = n_hashes;
        if ( summed_count >= threshold_cutoff )
        {
            // ensure that count was not already found for target with higher count
            // can happen in case of ambiguos targets in multiple filters
            if ( summed_count > std::get< 0 >( matches[target] ) )
            {
                matches[target] = std::make_tuple( summed_count, filter.filter_config.target_fpr[target] );
                if ( summed_count > max_count_read )
                    max_count_read = summed_count;
                if ( summed_count < min_count_read )
                    min_count_read = summed_count;
            }
        }
    }
}

void select_matches( Filter< THIBF >&       filter,
                     TMatches&              matches,
                     std::vector< size_t >& hashes,
                     auto&                  agent,
                     size_t                 threshold_cutoff,
                     size_t&                max_count_read,
                     size_t&                min_count_read,
                     size_t                 n_hashes )
{
    // Count only matches above threhsold
    seqan3::counting_vector< detail::TIntCount > counts = agent.bulk_count( hashes, threshold_cutoff );

    // Only one bin per target
    for ( auto const& [target, bins] : filter.map )
    {
        if ( counts[bins[0]] > 0 )
        {
            // Sum counts among bins (split target (user bins) into several technical bins)
            size_t summed_count = counts[bins[0]];
            // summed_count can be higher than n_hashes for matches in several technical bins
            if ( summed_count > n_hashes )
                summed_count = n_hashes;
            // ensure that count was not already found for target with higher count
            // can happen in case of ambiguous targets in multiple filters
            if ( summed_count > std::get< 0 >( matches[target] ) )
            {
                matches[target] = std::make_tuple( summed_count, filter.filter_config.target_fpr[target] );
                if ( summed_count > max_count_read )
                    max_count_read = summed_count;
                if ( summed_count < min_count_read )
                    min_count_read = summed_count;
            }
        }
    }
}

size_t filter_matches(
    ReadOut& read_out, TMatches& matches, TRep& rep, size_t n_hashes, double threshold_filter, double min_fpr_query )
{

    for ( auto const& [target, count_fpr] : matches )
    {
        if ( std::get< 0 >( count_fpr ) >= threshold_filter )
        {
            // Filter by fpr-query
            if ( min_fpr_query < 1.0 )
            {
                double q = 1;
                for ( size_t i = 0; i <= std::get< 0 >( count_fpr ); i++ )
                {
                    q -= binom( n_hashes, i ) * pow( std::get< 1 >( count_fpr ), i )
                         * pow( 1 - std::get< 1 >( count_fpr ), n_hashes - i );
                }
                if ( q > min_fpr_query )
                {
                    rep[{ read_out.prefix, target }].discarded_matches_fprquery++;
                    continue;
                }
            }

            rep[{ read_out.prefix, target }].matches++;
            read_out.matches.push_back( ReadMatch{ target, std::get< 0 >( count_fpr ) } );
        }
        else
        {
            rep[{ read_out.prefix, target }].discarded_matches_filter++;
        }
    }

    return read_out.matches.size();
}

void lca_matches( ReadOut& read_out, ReadOut& read_out_lca, size_t max_count_read, LCA& lca, TRep& rep )
{

    std::vector< std::string > targets;
    for ( auto const& r : read_out.matches )
    {
        targets.push_back( r.target );
    }

    std::string target_lca = lca.getLCA( targets );
    rep[{ read_out.prefix, target_lca }].seqs_lca++;
    read_out_lca.matches.push_back( ReadMatch{ target_lca, max_count_read } );
}


template < typename TFilter >
void classify( std::vector< Filter< TFilter > >& filters,
               LCA&                              lca,
               TRep&                             rep,
               TTotal&                           total,
               SafeQueue< ReadOut >&             classified_all_queue,
               SafeQueue< ReadOut >&             classified_lca_queue,
               SafeQueue< ReadOut >&             unclassified_queue,
               Config const&                     config,
               SafeQueue< ReadBatches >*         pointer_current,
               SafeQueue< ReadBatches >*         pointer_helper,
               HierarchyConfig const&            hierarchy_config,
               bool                              hierarchy_first,
               bool                              hierarchy_last )
{

    // oner hash adaptor per thread
    const auto minimiser_hash =
        seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ hierarchy_config.kmer_size } },
                                       seqan3::window_size{ hierarchy_config.window_size },
                                       seqan3::seed{ raptor::adjust_seed( hierarchy_config.kmer_size ) } );

    // one agent per thread per filter
    using TAgent = std::conditional_t< std::same_as< TFilter, THIBF >,
                                       THIBF::counting_agent_type< detail::TIntCount >,
                                       TIBF::counting_agent_type< detail::TIntCount > >;
    std::vector< TAgent > agents;
    for ( Filter< TFilter >& filter : filters )
    {
        agents.push_back( filter.ibf.template counting_agent< detail::TIntCount >() );
    }

    while ( true )
    {
        // Wait here until reads are available or push is over and queue is empty
        ReadBatches rb = pointer_current->pop();

        // If batch is empty exit thread
        if ( rb.ids.empty() )
            break;

        // store unclassified reads for next iteration
        ReadBatches left_over_reads{ rb.paired, rb.prefix };

        const size_t hashes_limit = std::numeric_limits< detail::TIntCount >::max();

        for ( size_t readID = 0; readID < rb.ids.size(); ++readID )
        {
            // read lenghts
            const size_t read1_len = rb.seqs[readID].size();
            const size_t read2_len = rb.paired ? rb.seqs2[readID].size() : 0;

            // Store matches for this read
            TMatches matches;

            // Best scoring kmer count
            size_t max_count_read = 0;
            size_t min_count_read = 0;
            size_t n_hashes       = 0;
            // if length is smaller than window, skip read
            if ( read1_len >= hierarchy_config.window_size )
            {
                // Count hashes
                std::vector< size_t > hashes = rb.seqs[readID] | minimiser_hash | seqan3::ranges::to< std::vector >();
                // Count hashes from both pairs if second is given
                if ( read2_len >= hierarchy_config.window_size )
                {
                    // Add hashes of second pair
                    const auto h2 = rb.seqs2[readID] | minimiser_hash | std::views::common;
                    hashes.insert( hashes.end(), h2.begin(), h2.end() );
                }

                n_hashes = hashes.size();
                // set min as max. possible hashes
                min_count_read = n_hashes;
                // if n_hashes are bigger than int limit, skip read
                if ( n_hashes <= hashes_limit )
                {
                    // Sum sequence to totals
                    if ( hierarchy_first )
                    {
                        total[rb.prefix].seqs_processed++;
                        total[rb.prefix].length_processed += read1_len + read2_len;
                        total[rb.prefix].kmers_processed += n_hashes;
                    }

                    // For each filter in the hierarchy
                    for ( size_t i = 0; i < filters.size(); ++i )
                    {
                        // Calculate threshold for cutoff (keep matches above)
                        size_t threshold_cutoff = threshold_rel( n_hashes, filters[i].filter_config.rel_cutoff );

                        // reset low threshold_cutoff to just one kmer (0 would match everywhere)
                        if ( threshold_cutoff == 0 )
                            threshold_cutoff = 1;

                        // count and select matches
                        select_matches( filters[i],
                                        matches,
                                        hashes,
                                        agents[i],
                                        threshold_cutoff,
                                        max_count_read,
                                        min_count_read,
                                        n_hashes );
                    }
                }
                else
                {
                    if ( hierarchy_first )
                        total[rb.prefix].seqs_skipped_big++;
                }
            }
            else
            {
                if ( hierarchy_first )
                    total[rb.prefix].seqs_skipped_small++;
            }

            // store read and matches to be printed
            ReadOut read_out( rb.ids[readID], rb.prefix );

            // if read got valid matches (above cutoff)
            if ( max_count_read > 0 )
            {

                // Calculate threshold for filtering (keep matches above)
                const size_t threshold_filter =
                    max_count_read - threshold_rel( max_count_read - min_count_read, hierarchy_config.rel_filter );

                // Filter matches
                const size_t count_filtered_matches =
                    filter_matches( read_out, matches, rep, n_hashes, threshold_filter, hierarchy_config.fpr_query );

                if ( count_filtered_matches > 0 )
                {
                    total[rb.prefix].seqs_classified++;
                    total[rb.prefix].kmers_from_classified_seqs += n_hashes;
                    total[rb.prefix].kmers_matches += max_count_read;

                    if ( !config.skip_lca )
                    {
                        ReadOut read_out_lca( rb.ids[readID], rb.prefix );
                        if ( count_filtered_matches == 1 )
                        {
                            // just one match, copy read read_out and set as unique
                            read_out_lca = read_out;
                            rep[{ read_out.prefix, read_out.matches[0].target }].seqs_unique++;
                        }
                        else
                        {
                            lca_matches( read_out, read_out_lca, max_count_read, lca, rep );
                        }

                        if ( config.output_lca )
                            classified_lca_queue.push( read_out_lca );
                    }
                    else
                    {
                        // Not running lca and has unique match
                        if ( count_filtered_matches == 1 )
                        {
                            rep[{ read_out.prefix, read_out.matches[0].target }].seqs_unique++;
                        }
                        else
                        {
                            // without tax, no lca, count multi-matches to a root node
                            // to keep consistency among reports (no. of classified reads)
                            rep[{ read_out.prefix, config.tax_root_node }].seqs_lca++;
                        }
                    }

                    if ( config.output_all )
                        classified_all_queue.push( read_out );

                    // read classified, continue to the next
                    continue;
                }
            }

            // not classified
            if ( !hierarchy_last ) // if there is more levels, store read
            {
                left_over_reads.ids.push_back( std::move( rb.ids[readID] ) );
                left_over_reads.seqs.push_back( std::move( rb.seqs[readID] ) );

                if ( rb.paired )
                {
                    left_over_reads.seqs2.push_back( std::move( rb.seqs2[readID] ) );
                }
            }
            else if ( config.output_unclassified ) // no more levels and no classification, add to
                                                   // unclassified printing queue
            {
                unclassified_queue.push( read_out );
            }
        }

        // if there are more levels to classify and something was left, keep reads in memory
        if ( !hierarchy_last && left_over_reads.ids.size() > 0 )
            pointer_helper->push( std::move( left_over_reads ) );
    }
}

void write_report( TRep& rep, TTax& tax, std::map< std::string, std::ofstream >& out, std::string hierarchy_label )
{
    for ( auto const& [idx, report] : rep )
    {
        auto prefix = idx.first;
        auto target = idx.second;

        if ( report.matches || report.seqs_lca || report.seqs_unique )
        {
            out[prefix] << hierarchy_label << '\t' << target << '\t' << report.matches << '\t' << report.seqs_unique
                        << '\t' << report.seqs_lca;

            if ( !tax.empty() )
            {
                out[prefix] << '\t' << tax.at( target ).rank << '\t' << tax.at( target ).name;
            }
            out[prefix] << '\n';
        }
    }
}

void write_report_totals( Stats& stats, std::map< std::string, std::ofstream >& out )
{
    for ( auto const& [prefix, total] : stats.total )
    {
        out[prefix] << "#total_classified\t" << total.seqs_classified << '\n';
        // account for unclassified and skipped sequences
        out[prefix] << "#total_unclassified\t" << total.input_seqs - total.seqs_classified << '\n';
    }
}

static inline void replace_all( std::string& str, const std::string& from, const std::string& to )
{
    size_t start_pos = 0;
    while ( ( start_pos = str.find( from, start_pos ) ) != std::string::npos )
    {
        str.replace( start_pos, from.length(), to );
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
}

size_t load_filter( THIBF&             filter,
                    IBFConfig&         ibf_config,
                    TBinMap&           bin_map,
                    std::string const& input_filter_file,
                    TTargetFpr&        target_fpr )
{
    std::ifstream              is( input_filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );

    uint32_t                                  parsed_version;
    uint64_t                                  window_size;
    seqan3::shape                             shape{};
    uint8_t                                   parts;
    bool                                      compressed;
    std::vector< std::vector< std::string > > bin_path{};
    double                                    fpr{};
    bool                                      is_hibf{};

    archive( parsed_version );
    archive( window_size );
    archive( shape );
    archive( parts );
    archive( compressed );
    archive( bin_path );
    archive( fpr );
    archive( is_hibf );
    archive( filter );

    // load ibf_config from raptor params
    ibf_config.window_size = window_size;
    ibf_config.kmer_size   = shape.count();
    ibf_config.max_fp      = fpr;

    // Create map from paths
    size_t binno{};
    for ( auto const& file_list : bin_path )
    {
        for ( auto const& filename : file_list )
        {
            // based on the filename, get target.minimiser
            // (e.g. 562.minimiser or GCF_013391805.1.minimiser), otherwise use filename as target
            auto   f     = std::filesystem::path( filename ).filename().string();
            size_t found = f.find( ".minimiser" );
            if ( found != std::string::npos )
                f = f.substr( 0, found );

            // workaround when file has a . (e.g. GCF_013391805.1)
            // "." replaced by "|||" in ganon build wrapper
            // fixed on ganon v2.0.0 (+ raptor 3.0.1) but kept for compatibility (>= ganon v1.8.0)
            replace_all( f, "|||", "." );

            // workaround when target has a space (e.g. s__Pectobacterium carotovorum)
            // " " replaced by "---" in ganon build wrapper
            replace_all( f, "---", " " );

            bin_map.push_back( std::make_tuple( binno, f ) );
            // same fpr for all
            target_fpr[f] = fpr;
        }
        ++binno;
    }

    return filter.user_bins.num_user_bins();
}

inline double false_positive( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t n_hashes )
{
    /*
     * calculates the theoretical false positive of a bin (bf) based on parameters
     */
    return std::pow( 1 - std::exp( -hash_functions / ( bin_size_bits / static_cast< double >( n_hashes ) ) ),
                     hash_functions );
}

size_t load_filter( TIBF&              filter,
                    IBFConfig&         ibf_config,
                    TBinMap&           bin_map,
                    std::string const& input_filter_file,
                    TTargetFpr&        target_fpr )
{
    std::ifstream              is( input_filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );

    std::tuple< int, int, int >                        parsed_version;
    std::vector< std::tuple< std::string, uint64_t > > hashes_count_std;

    archive( parsed_version );
    archive( ibf_config );
    archive( hashes_count_std );
    archive( bin_map );
    archive( filter );


    // generate fpr for each bin
    for ( auto const& [target, count] : hashes_count_std )
    {
        // Use average number of hashes for each bin to calculate fp
        uint64_t n_bins_target = std::ceil( count / static_cast< double >( ibf_config.max_hashes_bin ) );
        // this can be off by a very small number (rounding ceil on multiple bins)
        uint64_t n_hashes_bin = std::ceil( count / static_cast< double >( n_bins_target ) );

        // false positive for the current target, considering split bins
        target_fpr[target] =
            1.0
            - std::pow( 1.0 - false_positive( ibf_config.bin_size_bits, ibf_config.hash_functions, n_hashes_bin ),
                        n_bins_target );
        ;
    }


    return filter.bin_count();
}

TTax load_tax( std::string tax_file )
{
    TTax          tax;
    std::string   line;
    std::ifstream infile;
    infile.open( tax_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field );
        tax[fields[0]] = Node{ fields[1], fields[2], fields[3] };
    }
    infile.close();
    return tax;
}

template < typename TFilter >
bool load_files( std::vector< Filter< TFilter > >& filters, std::vector< FilterConfig >& fconf )
{
    size_t filter_cnt = 0;
    for ( auto& filter_config : fconf )
    {
        TTax       tax;
        IBFConfig  ibf_config;
        TBinMap    bin_map;
        TFilter    filter;
        TTargetFpr target_fpr;
        auto       bin_count = load_filter( filter, ibf_config, bin_map, filter_config.ibf_file, target_fpr );

        // Parse vector with bin_map to the old map
        TMap map;
        for ( auto const& [binno, target] : bin_map )
        {
            map[target].push_back( binno );
        }

        filter_config.ibf_config = ibf_config;
        filter_config.target_fpr = target_fpr;

        if ( filter_config.tax_file != "" )
            tax = load_tax( filter_config.tax_file );

        filters.push_back(
            Filter< TFilter >{ std::move( filter ), std::move( map ), std::move( tax ), bin_count, filter_config } );
        filter_cnt++;
    }

    return true;
}

void print_time( const StopClock& timeGanon, const StopClock& timeLoadFilters, const StopClock& timeClassPrint )
{
    using ::operator<<;
    std::cerr << "ganon-classify        start time: " << StopClock_datetime( timeGanon.begin() ) << std::endl;
    std::cerr << "ganon-classify          end time: " << StopClock_datetime( timeGanon.end() ) << std::endl;
    std::cerr << "loading filter(s)    elapsed (s): " << timeLoadFilters.elapsed() << " seconds" << std::endl;
    std::cerr << "classifying+printing elapsed (s): " << timeClassPrint.elapsed() << " seconds" << std::endl;
    std::cerr << "total                elapsed (s): " << timeGanon.elapsed() << " seconds" << std::endl;
    std::cerr << "----------------------------------------------------------------------" << std::endl;
    std::cerr << std::endl;
}

void print_stats_db( const Total& total, double seq_processed, size_t seq_unclassified )
{
    const size_t seq_multiple_matches = total.seqs_classified - total.seqs_unique;
    const double avg_seq_matches =
        total.seqs_classified ? ( total.matches / static_cast< double >( total.seqs_classified ) ) : 0;
    const double kmers_matched_perc =
        total.kmers_matches ? ( total.kmers_matches / static_cast< double >( total.kmers_from_classified_seqs ) ) * 100
                            : 0;

    // std::cerr << std::fixed << std::setprecision( 4 );
    std::cerr << "" << total.seqs_classified << " sequences classified ("
              << ( total.seqs_classified / seq_processed ) * 100 << "%)" << std::endl;
    std::cerr << "  " << total.seqs_unique << " with unique matches (" << ( total.seqs_unique / seq_processed ) * 100
              << "%)" << std::endl;
    std::cerr << "  " << seq_multiple_matches << " with multiple matches ("
              << ( seq_multiple_matches / seq_processed ) * 100 << "%)" << std::endl;
    if ( seq_unclassified > 0 )
    {
        std::cerr << "" << seq_unclassified << " sequences unclassified (" << ( seq_unclassified / seq_processed ) * 100
                  << "%)" << std::endl;
        if ( total.seqs_skipped_small )
        {
            std::cerr << "  " << total.seqs_skipped_small << " sequences skipped (shorter than window size)"
                      << std::endl;
        }
        if ( total.seqs_skipped_big )
        {
            std::cerr << "  " << total.seqs_skipped_big
                      << " sequences skipped (larger than allowed, check compilation with -DLONGREADS)" << std::endl;
        }
    }
    std::cerr << "matches: " << total.matches << " (avg. " << avg_seq_matches << " reference/sequence), "
              << total.discarded_matches_filter << " discarded (--rel-filter), " << total.discarded_matches_fprquery
              << " discarded (--fpr-query)" << std::endl;
    std::cerr << "k-mers: " << total.kmers_matches << "/" << total.kmers_from_classified_seqs
              << " k-mers matched/k-mers from classified sequences" << " (" << kmers_matched_perc << "%)" << std::endl;
}

void print_stats( Stats&                                          stats,
                  double                                          elapsed_classification,
                  const std::map< std::string, HierarchyConfig >& parsed_hierarchy )
{
    std::cerr << "ganon-classify processed " << stats.total_seqs_processed << " sequences ("
              << stats.total_length_processed / 1000000.0 << " Mbp) with " << stats.total_kmers_processed
              << " k-mers in " << elapsed_classification << " seconds ("
              << ( stats.total_length_processed / 1000000.0 ) / ( elapsed_classification / 60.0 ) << " Mbp/m)"
              << std::endl;

    for ( auto const& [prefix, total] : stats.total )
    {
        if ( stats.total.size() > 1 )
        {
            std::cerr << std::endl;
            std::cerr << "[" << prefix << "] " << total.seqs_processed << " sequences ("
                      << total.length_processed / 1000000.0 << " Mbp) with " << total.kmers_processed << " k-mers"
                      << std::endl;
        }

        const size_t seq_unclassified = total.seqs_processed - total.seqs_classified;
        const double seq_processed    = total.seqs_processed > 0 ? static_cast< double >( total.seqs_processed )
                                                                 : 1; // to not report nan on divisions
        detail::print_stats_db( total, seq_processed, seq_unclassified );

        if ( parsed_hierarchy.size() > 1 )
        {
            std::cerr << std::endl;
            std::cerr << "By database hierarchical level:" << std::endl;
            for ( auto const& h : parsed_hierarchy )
            {
                std::string hierarchy_label = h.first;
                std::cerr << hierarchy_label << ":" << std::endl;
                detail::print_stats_db( stats.hierarchy_total[hierarchy_label][prefix], seq_processed, 0 );
            }
        }
    }
}

void write_stats_db( const Total&   total,
                     double         seq_processed,
                     size_t         seq_unclassified,
                     size_t         kmers_processed,
                     std::string    prefix,
                     std::string    hierarchy_level,
                     std::ofstream& out_sta )
{

    const size_t seq_multiple_matches = total.seqs_classified - total.seqs_unique;
    const double avg_seq_matches =
        total.seqs_classified ? ( total.matches / static_cast< double >( total.seqs_classified ) ) : 0;
    const double kmers_matched_perc =
        total.kmers_matches ? ( total.kmers_matches / static_cast< double >( total.kmers_from_classified_seqs ) ) * 100
                            : 0;

    out_sta << std::fixed << std::setprecision( 6 );
    out_sta << prefix << '\t';
    out_sta << hierarchy_level << '\t';
    out_sta << static_cast< size_t >( seq_processed ) << '\t';
    out_sta << seq_unclassified << '\t';
    out_sta << total.seqs_classified << '\t';
    out_sta << ( total.seqs_classified / seq_processed ) * 100 << '\t';
    out_sta << total.seqs_unique << '\t';
    out_sta << ( total.seqs_unique / seq_processed ) * 100 << '\t';
    out_sta << seq_multiple_matches << '\t';
    out_sta << ( seq_multiple_matches / seq_processed ) * 100 << '\t';
    out_sta << total.matches << '\t';
    out_sta << avg_seq_matches << '\t';
    out_sta << total.discarded_matches_filter << '\t';
    out_sta << total.discarded_matches_fprquery << '\t';
    out_sta << kmers_processed << '\t';
    out_sta << total.kmers_matches << '\t';
    out_sta << total.kmers_from_classified_seqs << '\t';
    out_sta << kmers_matched_perc << '\n';
}

void write_stats( std::string                                     output_prefix,
                  Stats&                                          stats,
                  const std::map< std::string, HierarchyConfig >& parsed_hierarchy )
{

    for ( auto const& [prefix, total] : stats.total )
    {
        std::ofstream out_sta{ output_prefix + prefix + ".sta" };

        out_sta << "prefix" << '\t';
        out_sta << "hierarchy_label" << '\t';
        out_sta << "seq_processed" << '\t';
        out_sta << "seq_unclassified" << '\t';
        out_sta << "seq_classified" << '\t';
        out_sta << "seq_classified_perc" << '\t';
        out_sta << "seq_unique_matches" << '\t';
        out_sta << "seq_unique_matches_perc" << '\t';
        out_sta << "seq_multiple_matches" << '\t';
        out_sta << "seq_multiple_matches_perc" << '\t';
        out_sta << "matches" << '\t';
        out_sta << "avg_matches_ref_seq" << '\t';
        out_sta << "dis_matches_rel_filter" << '\t';
        out_sta << "dis_matches_fpr_query" << '\t';
        out_sta << "kmers_proccessed" << '\t';
        out_sta << "kmers_matched" << '\t';
        out_sta << "kmers_from_classified_seqs" << '\t';
        out_sta << "kmers_matched_perc" << '\n';

        const size_t seq_unclassified = total.seqs_processed - total.seqs_classified;
        const double seq_processed    = total.seqs_processed > 0 ? static_cast< double >( total.seqs_processed )
                                                                 : 1; // to not report nan on divisions

        for ( auto const& h : parsed_hierarchy )
        {
            detail::write_stats_db( stats.hierarchy_total[h.first][prefix],
                                    seq_processed,
                                    seq_unclassified,
                                    total.kmers_processed,
                                    prefix,
                                    h.first,
                                    out_sta );
        }

        // If multiple hierarchy levels, print -total- line with sum
        if ( parsed_hierarchy.size() > 1 )
        {
            detail::write_stats_db(
                total, seq_processed, seq_unclassified, total.kmers_processed, prefix, "-total-", out_sta );
        }
        out_sta.close();
    }
}

void parse_reads( SafeQueue< ReadBatches >& queue1, Stats& stats, Config const& config, TReadConfig& reads_config )
{

    for ( auto& [prefix, files] : reads_config )
    {
        for ( auto& [filename1, filename2] : files )
        {
            bool paired = filename2.empty() ? false : true;
            try
            {
                if ( paired )
                {
                    seqan3::sequence_file_input< raptor::dna4_traits,
                                                 seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                        fin1{ filename1 };

                    seqan3::sequence_file_input< raptor::dna4_traits,
                                                 seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                        fin2{ filename2 };

                    for ( auto&& rec : fin1 | seqan3::views::chunk( config.n_reads ) )
                    {
                        ReadBatches rb{ true };
                        for ( auto& [id, seq] : rec )
                        {
                            rb.ids.push_back( std::move( id ) );
                            rb.seqs.push_back( std::move( seq ) );
                        }
                        // loop in the second file and get same amount of reads
                        for ( auto& [id, seq] : fin2 | std::views::take( config.n_reads ) )
                        {
                            rb.seqs2.push_back( std::move( seq ) );
                        }
                        stats.total[prefix].input_seqs += rb.ids.size();
                        rb.prefix = prefix;
                        queue1.push( std::move( rb ) );
                    }
                }
                else
                {
                    seqan3::sequence_file_input< raptor::dna4_traits,
                                                 seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                        fin1{ filename1 };

                    for ( auto&& rec : fin1 | seqan3::views::chunk( config.n_reads ) )
                    {
                        ReadBatches rb{ false };
                        for ( auto& [id, seq] : rec )
                        {
                            rb.ids.push_back( std::move( id ) );
                            rb.seqs.push_back( std::move( seq ) );
                        }
                        stats.total[prefix].input_seqs += rb.ids.size();
                        rb.prefix = prefix;
                        queue1.push( std::move( rb ) );
                    }
                }
            }
            catch ( seqan3::parse_error const& ext )
            {
                std::cerr << "Error parsing file(s) [" << filename1 << ", " << filename2 << "]" << ext.what()
                          << std::endl;
                continue;
            }
        }
    }
    queue1.notify_push_over();
}

void write_classified( SafeQueue< ReadOut >& classified_queue, std::map< std::string, std::ofstream >& out )
{
    while ( true )
    {
        ReadOut ro = classified_queue.pop();
        if ( ro.readID != "" )
        {
            for ( size_t i = 0; i < ro.matches.size(); ++i )
            {
                out[ro.prefix] << ro.readID << '\t' << ro.matches[i].target << '\t' << ro.matches[i].kmer_count << '\n';
            }
        }
        else
        {
            break;
        }
    }
}

void write_unclassified( SafeQueue< ReadOut >& unclassified_queue, std::map< std::string, std::ofstream >& out )
{
    while ( true )
    {
        ReadOut ro = unclassified_queue.pop();
        if ( ro.readID != "" )
        {
            out[ro.prefix] << ro.readID << '\n';
        }
        else
        {
            break;
        }
    }
}

template < typename TFilter >
TTax merge_tax( std::vector< Filter< TFilter > > const& filters )
{
    if ( filters.size() == 1 )
    {
        return filters[0].tax;
    }
    else
    {
        TTax merged_tax = filters[0].tax;
        for ( size_t i = 1; i < filters.size(); ++i )
        {
            // merge taxonomies keeping the first one as a default
            merged_tax.insert( filters[i].tax.begin(), filters[i].tax.end() );
        }
        return merged_tax;
    }
}

template < typename TFilter >
void validate_targets_tax( std::vector< Filter< TFilter > > const& filters,
                           TTax&                                   tax,
                           bool                                    quiet,
                           const std::string                       tax_root_node )
{
    for ( auto const& filter : filters )
    {
        for ( auto const& [target, bins] : filter.map )
        {
            if ( tax.count( target ) == 0 )
            {
                tax[target] = Node{ tax_root_node, "no rank", target };
                if ( !quiet )
                    std::cerr << "WARNING: target [" << target << "] without tax entry, setting parent as root node ["
                              << tax_root_node << "]" << std::endl;
            }
        }
    }
}

void pre_process_lca( LCA& lca, TTax& tax, std::string tax_root_node )
{
    for ( auto const& [target, node] : tax )
    {
        lca.addEdge( node.parent, target );
    }
    lca.doEulerWalk( tax_root_node );
}

} // namespace detail

template < typename TFilter >
bool ganon_classify( Config config )
{
    // Start time count
    StopClock timeGanon;
    timeGanon.start();

    auto parsed_hierarchy = detail::parse_hierarchy( config );

    // Prepare reads
    detail::TReadConfig reads_config;
    if ( !detail::parse_reads_config( config, reads_config ) )
        return false;

    // Create dirs if necessary
    for ( auto& [prefix, files] : reads_config )
    {
        std::filesystem::path filepath = std::string( config.output_prefix + prefix );
        if ( !std::filesystem::is_directory( filepath ) && !filepath.parent_path().empty() )
        {
            std::filesystem::create_directories( filepath.parent_path() );
        }
    }

    if ( config.verbose )
    {
        detail::print_hierarchy( parsed_hierarchy );
        detail::print_reads_config( reads_config );
        detail::print_output_files( config, parsed_hierarchy, reads_config );
    }

    // Initialize variables
    StopClock timeLoadFilters;
    StopClock timeClassPrint;

    detail::Stats                          stats;
    std::map< std::string, std::ofstream > out_rep; // Set default output stream (file or stdout)
    std::map< std::string, std::ofstream > out_all; // output all file by read prefix
    std::map< std::string, std::ofstream > out_lca; // output lca file
    std::map< std::string, std::ofstream > out_unc; // output unclassified file

    for ( auto& [prefix, files] : reads_config )
    {
        out_rep[prefix].open( config.output_prefix + prefix + ".rep" );
    }

    // Queues for internal read handling
    // queue1 get reads from file
    // queue2 will get unclassified reads if hierachy == 2
    // if hierachy == 3 queue1 is used for unclassified and so on
    // config.n_batches*config.n_reads = max. amount of reads in memory
    SafeQueue< detail::ReadBatches >  queue1( config.n_batches );
    SafeQueue< detail::ReadBatches >  queue2;
    SafeQueue< detail::ReadBatches >* pointer_current = &queue1; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_helper  = &queue2; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_extra;             // pointer to the queues

    // Define one threads for decompress bgzf files
    seqan3::contrib::bgzf_thread_count = 1u;

    // Thread for reading input files
    std::future< void > read_task = std::async( std::launch::async,
                                                detail::parse_reads,
                                                std::ref( queue1 ),
                                                std::ref( stats ),
                                                std::ref( config ),
                                                std::ref( reads_config ) );

    // Thread for printing unclassified reads
    SafeQueue< detail::ReadOut > unclassified_queue;
    std::future< void >          write_unclassified_task;
    if ( config.output_unclassified )
    {
        for ( auto& [prefix, files] : reads_config )
        {
            out_unc[prefix].open( config.output_prefix + prefix + ".unc" );
        }

        write_unclassified_task = std::async(
            std::launch::async, detail::write_unclassified, std::ref( unclassified_queue ), std::ref( out_unc ) );
    }


    // Classify reads iteractively for each hierarchy level
    size_t       hierarchy_id   = 0;
    const size_t hierarchy_size = parsed_hierarchy.size();
    for ( auto& [hierarchy_label, hierarchy_config] : parsed_hierarchy )
    {
        ++hierarchy_id;
        bool                                     hierarchy_first = ( hierarchy_id == 1 );
        bool                                     hierarchy_last  = ( hierarchy_id == hierarchy_size );
        std::vector< detail::Filter< TFilter > > filters;
        detail::TTax                             tax;
        LCA                                      lca;

        timeLoadFilters.start();
        const bool loaded = detail::load_files( filters, parsed_hierarchy[hierarchy_label].filters );
        if ( !loaded )
        {
            std::cerr << "ERROR: loading ibf or tax files" << std::endl;
            return false;
        }
        timeLoadFilters.stop();

        hierarchy_config.kmer_size   = filters[0].filter_config.ibf_config.kmer_size;
        hierarchy_config.window_size = filters[0].filter_config.ibf_config.window_size;
        if ( filters.size() > 1 )
        {
            // Check if all filters share the same k,w
            for ( auto const& filter : filters )
            {
                if ( filter.filter_config.ibf_config.kmer_size != hierarchy_config.kmer_size
                     || filter.filter_config.ibf_config.window_size != hierarchy_config.window_size )
                {
                    std::cerr << "ERROR: databases on the same hierarchy should share same k-mer and window sizes"
                              << std::endl;
                    return false;
                }
            }
        }


        // Parse tax if provided
        if ( filters[0].filter_config.tax_file != "" )
        {
            // merge repeated elements from multiple filters
            tax = detail::merge_tax( filters );
            // if target not found in tax, add node target with parent root
            detail::validate_targets_tax( filters, tax, config.quiet, config.tax_root_node );
        }

        if ( !config.skip_lca )
        {
            if ( tax.count( config.tax_root_node ) == 0 )
            {
                std::cerr << "Root node [" << config.tax_root_node << "] not found (--tax-root-node)" << std::endl;
                return false;
            }
            // pre-processing of nodes to generate LCA
            detail::pre_process_lca( lca, tax, config.tax_root_node );
        }


        // Thread for printing classified reads (.one, .all)
        std::vector< std::future< void > > write_tasks;

        // hierarchy_id = 1
        //  pointer_current=queue1, data comes from file in a limited size queue
        //  pointer_helper=queue2, empty
        // hierarchy_id > 1
        //  pointer_current=queue2, with all data already in from last iteration
        //  pointer_helper=queue1, empty
        // Exchange queues instance pointers for each hierachy (if not first)
        if ( !hierarchy_first )
        {
            pointer_extra   = pointer_current;
            pointer_current = pointer_helper;
            pointer_helper  = pointer_extra;

            // Remove size limit from reading since it's always already loaded
            if ( hierarchy_id == 2 )
                queue1.set_max_size( std::numeric_limits< size_t >::max() );
        }

        SafeQueue< detail::ReadOut > classified_all_queue;
        SafeQueue< detail::ReadOut > classified_lca_queue;
        // Append to files if not first hierarchy or not output single
        const auto file_mode = hierarchy_first || !config.output_single ? std::ofstream::out : std::ofstream::app;

        if ( config.output_lca && !config.skip_lca )
        {
            for ( auto& [prefix, files] : reads_config )
            {
                // open file once for each read prefix
                out_lca[prefix].open( config.output_prefix + prefix + "." + hierarchy_config.output_file_lca,
                                      file_mode );
            }

            // Start writing thread for lca matches
            write_tasks.emplace_back( std::async(
                std::launch::async, detail::write_classified, std::ref( classified_lca_queue ), std::ref( out_lca ) ) );
        }
        if ( config.output_all )
        {
            for ( auto& [prefix, files] : reads_config )
            {
                // open file once for each read prefix, append if output_single and not first
                out_all[prefix].open( config.output_prefix + prefix + "." + hierarchy_config.output_file_all,
                                      file_mode );
            }

            // Start writing thread for all matches
            write_tasks.emplace_back( std::async(
                std::launch::async, detail::write_classified, std::ref( classified_all_queue ), std::ref( out_all ) ) );
        }


        // One report and total counters for each thread
        std::vector< detail::TRep >   reports( config.threads );
        std::vector< detail::TTotal > totals( config.threads );

        std::vector< std::future< void > > tasks;
        // Threads for classification
        timeClassPrint.start();
        for ( size_t taskNo = 0; taskNo < config.threads; ++taskNo )
        {

            tasks.emplace_back( std::async( std::launch::async,
                                            detail::classify< TFilter >,
                                            std::ref( filters ),
                                            std::ref( lca ),
                                            std::ref( reports[taskNo] ),
                                            std::ref( totals[taskNo] ),
                                            std::ref( classified_all_queue ),
                                            std::ref( classified_lca_queue ),
                                            std::ref( unclassified_queue ),
                                            std::ref( config ),
                                            pointer_current,
                                            pointer_helper,
                                            hierarchy_config,
                                            hierarchy_first,
                                            hierarchy_last ) );
        }

        // Wait here until classification is over
        for ( auto&& task : tasks )
        {
            task.get();
        }

        // After classification, no more reads are going to be pushed to the output
        classified_all_queue.notify_push_over();
        classified_lca_queue.notify_push_over();

        // Sum reports of each threads into one
        detail::TRep rep = sum_reports( reports );

        // Add reports values into stats
        stats.add_totals( hierarchy_label, totals );
        stats.add_reports( hierarchy_label, rep );

        // write reports
        detail::write_report( rep, tax, out_rep, hierarchy_label );

        // Wait here until all files are written
        for ( auto&& task : write_tasks )
        {
            task.get();
        }
        timeClassPrint.stop();

        if ( config.output_lca )
            for ( auto& [prefix, file] : out_lca )
                file.close();
        if ( config.output_all )
            for ( auto& [prefix, file] : out_all )
                file.close();

        if ( hierarchy_first )
        {
            read_task.get();                    // get reading tasks at the end of the first hierarchy
            pointer_helper->notify_push_over(); // notify push is over, only on first time (will be always set over for
                                                // next iterations)
        }
    }

    // Wait here until all unclassified reads are written
    if ( config.output_unclassified )
    {
        unclassified_queue.notify_push_over();
        write_unclassified_task.get();
    }
    if ( config.output_unclassified )
        for ( auto& [prefix, file] : out_unc )
            file.close();

    detail::write_report_totals( stats, out_rep );

    for ( auto& [prefix, file] : out_rep )
        file.close();

    if ( config.output_stats )
    {
        detail::write_stats( config.output_prefix, stats, parsed_hierarchy );
    }

    timeGanon.stop();

    if ( !config.quiet )
    {
        if ( config.verbose )
        {
            detail::print_time( timeGanon, timeLoadFilters, timeClassPrint );
        }

        detail::print_stats( stats, timeClassPrint.elapsed(), parsed_hierarchy );
    }

    return true;
}

bool run( Config config )
{

    // Validate configuration input
    if ( !config.validate() )
        return false;

    // Print config
    if ( config.verbose )
        std::cerr << config;

    if ( config.hibf )
        return ganon_classify< detail::THIBF >( config );
    else
        return ganon_classify< detail::TIBF >( config );
}

} // namespace GanonClassify
