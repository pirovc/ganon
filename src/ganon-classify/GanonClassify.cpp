#include "GanonClassify.hpp"
#include "hierarchical_interleaved_bloom_filter.hpp"

#include <robin_hood.h>

#include <utils/LCA.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>
#include <utils/adjust_seed.hpp>
#include <utils/dna4_traits.hpp>

#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <cinttypes>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace GanonClassify
{

namespace detail
{

typedef raptor::hierarchical_interleaved_bloom_filter< seqan3::data_layout::uncompressed > THIBF;
typedef seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed >              TIBF;
typedef robin_hood::unordered_map< std::string, uint16_t >                                 TMatches;

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
};

struct ReadMatch
{
    std::string target;
    uint16_t    kmer_count;
};

struct ReadOut
{
    ReadOut()
    {
    }

    ReadOut( std::string _readID )
    {
        readID = _readID;
    }

    std::string              readID;
    std::vector< ReadMatch > matches;
};

struct Rep
{
    // Report with counts of matches and reads assigned (unique or lca) for each target
    uint64_t matches      = 0;
    uint64_t lca_reads    = 0;
    uint64_t unique_reads = 0;
};

typedef robin_hood::unordered_map< std::string, Rep >      TRep;
typedef robin_hood::unordered_map< std::string, Node >     TTax;
typedef robin_hood::unordered_map< uint32_t, std::string > TMap;

struct Total
{
    uint64_t reads_processed  = 0;
    uint64_t length_processed = 0;
    uint64_t reads_classified = 0;
    uint64_t matches          = 0;
    uint64_t unique_matches   = 0;
};

struct Stats
{
    Total total;
    // number of reads in the input files
    uint64_t input_reads = 0;
    // Total for each hierarchy
    std::map< std::string, Total > hierarchy_total;

    void add_totals( std::string hierarchy_label, std::vector< Total > const& totals )
    {
        // add several totals (from threads) into the stats
        for ( auto const& t : totals )
        {
            total.reads_processed += t.reads_processed;
            total.length_processed += t.length_processed;
            total.reads_classified += t.reads_classified;
            hierarchy_total[hierarchy_label].reads_processed += t.reads_processed;
            hierarchy_total[hierarchy_label].reads_classified += t.reads_classified;
            hierarchy_total[hierarchy_label].length_processed += t.reads_classified;
        }
    }

    void add_reports( std::string hierarchy_label, TRep const& report )
    {
        // add values from reports to stats
        for ( auto const& [target, rep] : report )
        {
            total.matches += rep.matches;
            total.unique_matches += rep.unique_reads;
            hierarchy_total[hierarchy_label].matches += rep.matches;
            hierarchy_total[hierarchy_label].unique_matches += rep.unique_reads;
        }
    }
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


inline TRep sum_reports( std::vector< TRep > const& reports )
{
    TRep report_sum;
    for ( auto const& report : reports )
    {
        for ( auto const& [target, r] : report )
        {
            report_sum[target].matches += r.matches;
            report_sum[target].lca_reads += r.lca_reads;
            report_sum[target].unique_reads += r.unique_reads;
        }
    }
    return report_sum;
}

inline uint16_t threshold_abs( uint16_t kmers, uint16_t e, uint8_t k, uint8_t o )
{
    return ( kmers * o > e * k ) ? std::ceil( kmers - ( e * k ) / static_cast< double >( o ) ) : 0u;
}

inline uint16_t threshold_rel( uint16_t kmers, double p )
{
    return std::ceil( kmers * p );
}

inline uint16_t get_abs_error( uint16_t kmers, uint8_t k, uint8_t o, uint8_t count )
{
    return std::ceil( ( o * ( -count + kmers ) ) / static_cast< double >( k ) );
}

void select_matches( Filter< TIBF >&        filter,
                     TMatches&              matches,
                     std::vector< size_t >& hashes_f,
                     std::vector< size_t >& hashes_r,
                     auto&                  agent,
                     uint16_t               threshold_cutoff,
                     uint16_t&              max_kmer_count_read )
{
    // count matches
    seqan3::counting_vector< uint16_t > counts_f = agent.bulk_count( hashes_f );
    seqan3::counting_vector< uint16_t > counts_r = agent.bulk_count( hashes_r );

    // for each bin
    // for ( uint32_t bin_n = 0; bin_n < filter.ibf.noOfBins; ++bin_n )
    // loop in map structure to avoid extra validations when map.size() < filter.ibf.noOfBins when ibf is updated and
    // sequences removed also avoid the error of spurius results from empty bins (bug reported)
    for ( auto const& [bin_n, target] : filter.map )
    {
        // if kmer count is higher than threshold_cutoff
        if ( counts_f[bin_n] >= threshold_cutoff || counts_r[bin_n] >= threshold_cutoff )
        {
            // get best matching strand
            uint16_t max_kmer_count_bin = std::max( counts_f[bin_n], counts_r[bin_n] );
            // keep only the best match target/read when same targets are split in several
            // bins
            if ( matches.count( target ) == 0 || max_kmer_count_bin > matches[target] )
            {
                // store match to target
                matches[target] = max_kmer_count_bin;
                if ( max_kmer_count_bin > max_kmer_count_read )
                    max_kmer_count_read = max_kmer_count_bin;
            }
        }
    }
}

void select_matches( Filter< THIBF >&       filter,
                     TMatches&              matches,
                     std::vector< size_t >& hashes_f,
                     std::vector< size_t >& hashes_r,
                     auto&                  agent,
                     uint16_t               threshold_cutoff,
                     uint16_t&              max_kmer_count_read )
{
    // count matches, return only above threhsold
    seqan3::counting_vector< uint16_t > counts_f = agent.bulk_count( hashes_f, threshold_cutoff );

    // on HIBF, returns user bins
    for ( auto const& [bin_n, target] : filter.map )
    {
        if ( counts_f[bin_n] > 0 )
        {
            // keep check on target in case of multiple filters with same target
            if ( matches.count( target ) == 0 || counts_f[bin_n] > matches[target] )
            {
                // store match to target
                matches[target] = counts_f[bin_n];
                if ( counts_f[bin_n] > max_kmer_count_read )
                    max_kmer_count_read = counts_f[bin_n];
            }
        }
    }
}

auto get_offset_hashes( auto& hashes, uint8_t offset )
{
    // return offset of hashes, always keep last one to cover full read (extra_pos)
    uint16_t n_kmers     = hashes.size() - 1;
    auto     nth_element = [offset, n_kmers]( auto&& tuple ) {
        return ( std::get< 0 >( tuple ) % offset == 0 || n_kmers == std::get< 0 >( tuple ) );
    };
    return seqan3::views::zip( std::views::iota( 0 ), hashes ) | std::views::filter( nth_element ) | std::views::values
           | seqan3::views::to< std::vector >;
}

void get_hashes( std::vector< seqan3::dna4 >& seq,
                 std::vector< size_t >&       hashes_f,
                 std::vector< size_t >&       hashes_r,
                 uint8_t                      window_size,
                 uint8_t                      offset,
                 auto&                        kmer_hash,
                 auto&                        minimiser_hash )
{

    if ( window_size > 0 )
    {
        // minimizers - no need to calculate reverse hashes
        hashes_f = seq | minimiser_hash | seqan3::views::to< std::vector >;
    }
    else
    {
        if ( offset > 1 )
        {
            // offset
            auto hf  = seq | kmer_hash;
            hashes_f = get_offset_hashes( hf, offset );
            auto hr  = seq | std::views::reverse | seqan3::views::complement | kmer_hash;
            hashes_r = get_offset_hashes( hr, offset );
        }
        else
        {
            // kmer
            hashes_f = seq | kmer_hash | seqan3::views::to< std::vector >;
            hashes_r =
                seq | std::views::reverse | seqan3::views::complement | kmer_hash | seqan3::views::to< std::vector >;
        }
    }
}

void get_hashes( std::vector< seqan3::dna4 >& seq,
                 std::vector< seqan3::dna4 >& seq2,
                 std::vector< size_t >&       hashes_f,
                 std::vector< size_t >&       hashes_r,
                 uint8_t                      window_size,
                 uint8_t                      offset,
                 auto&                        kmer_hash,
                 auto&                        minimiser_hash )
{

    // first pair
    get_hashes( seq, hashes_f, hashes_r, window_size, offset, kmer_hash, minimiser_hash );

    // hashes are inserted in the same vector for the second pair in the rev.comp.
    if ( window_size > 0 )
    {
        // minimizers - no need to calculate reverse hashes
        auto hf = seq2 | minimiser_hash | std::views::common;
        hashes_f.insert( hashes_f.end(), hf.begin(), hf.end() );
    }
    else
    {
        if ( offset > 1 )
        {
            // offset
            auto hf        = seq2 | std::views::reverse | seqan3::views::complement | kmer_hash;
            auto hashes_f2 = get_offset_hashes( hf, offset );
            hashes_f.insert( hashes_f.end(), hashes_f2.begin(), hashes_f2.end() );
            auto hr        = seq2 | kmer_hash;
            auto hashes_r2 = get_offset_hashes( hr, offset );
            hashes_r.insert( hashes_r.end(), hashes_r2.begin(), hashes_r2.end() );
        }
        else
        {
            // kmer
            auto hf = seq2 | std::views::reverse | seqan3::views::complement | kmer_hash | std::views::common;
            hashes_f.insert( hashes_f.end(), hf.begin(), hf.end() );
            auto hr = seq2 | kmer_hash | std::views::common;
            hashes_r.insert( hashes_r.end(), hr.begin(), hr.end() );
        }
    }
}

uint16_t get_threshold_filter( HierarchyConfig const& hierarchy_config,
                               uint16_t               kmers,
                               uint16_t               max_kmer_count_read,
                               uint8_t                kmer_size,
                               uint8_t                offset )
{
    uint16_t threshold_filter = 0;
    if ( hierarchy_config.rel_filter >= 0 )
    {
        threshold_filter = max_kmer_count_read - threshold_rel( max_kmer_count_read, hierarchy_config.rel_filter );
    }
    else if ( hierarchy_config.abs_filter >= 0 )
    {
        // get maximum possible number of errors of best match + abs_filter
        uint16_t max_error_threshold =
            get_abs_error( kmers, kmer_size, offset, max_kmer_count_read ) + hierarchy_config.abs_filter;
        // get min kmer count necesary to achieve the calculated number of errors
        threshold_filter = threshold_abs( kmers, max_error_threshold, kmer_size, offset );
    }
    return threshold_filter;
}

uint32_t filter_matches( ReadOut& read_out, TMatches& matches, TRep& rep, uint16_t threshold_filter )
{

    for ( auto const& [target, kmer_count] : matches )
    {
        if ( kmer_count >= threshold_filter )
        {
            rep[target].matches++;
            read_out.matches.push_back( ReadMatch{ target, kmer_count } );
        }
    }

    return read_out.matches.size();
}

void lca_matches( ReadOut& read_out, ReadOut& read_out_lca, uint16_t max_kmer_count_read, LCA& lca, TRep& rep )
{

    std::vector< std::string > targets;
    for ( auto const& r : read_out.matches )
    {
        targets.push_back( r.target );
    }

    std::string target_lca = lca.getLCA( targets );
    rep[target_lca].lca_reads++;
    read_out_lca.matches.push_back( ReadMatch{ target_lca, max_kmer_count_read } );
}


template < typename TFilter >
void classify( std::vector< Filter< TFilter > >& filters,
               LCA&                              lca,
               TRep&                             rep,
               Total&                            total,
               SafeQueue< ReadOut >&             classified_all_queue,
               SafeQueue< ReadOut >&             classified_lca_queue,
               SafeQueue< ReadOut >&             unclassified_queue,
               Config const&                     config,
               SafeQueue< ReadBatches >*         pointer_current,
               SafeQueue< ReadBatches >*         pointer_helper,
               HierarchyConfig const&            hierarchy_config,
               bool                              hierarchy_first,
               bool                              hierarchy_last,
               bool                              run_lca )
{

    // get max from kmer/window size
    uint16_t wk_size = ( hierarchy_config.window_size > 0 ) ? hierarchy_config.window_size : hierarchy_config.kmer_size;

    // oner hash adaptor per thread
    auto kmer_hash = seqan3::views::kmer_hash( seqan3::ungapped{ hierarchy_config.kmer_size } );
    auto minimiser_hash =
        seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ hierarchy_config.kmer_size } },
                                       seqan3::window_size{ wk_size },
                                       seqan3::seed{ raptor::adjust_seed( hierarchy_config.kmer_size ) } );

    // one agent per thread per filter
    using TAgent = std::conditional_t< std::same_as< TFilter, THIBF >,
                                       THIBF::counting_agent_type< uint16_t >,
                                       TIBF::counting_agent_type< uint16_t > >;
    std::vector< TAgent > agents;
    for ( Filter< TFilter >& filter : filters )
    {
        agents.push_back( filter.ibf.counting_agent() );
    }

    while ( true )
    {
        // Wait here until reads are available or push is over and queue is empty
        ReadBatches rb = pointer_current->pop();

        // If batch is empty exit thread
        if ( rb.ids.empty() )
            break;

        // store unclassified reads for next iteration
        ReadBatches left_over_reads{ rb.paired };

        for ( uint32_t readID = 0; readID < rb.ids.size(); ++readID )
        {
            // read lenghts
            uint16_t read1_len = rb.seqs[readID].size();
            uint16_t read2_len = rb.paired ? rb.seqs2[readID].size() : 0;

            // Store matches for this read
            TMatches matches;

            // Retrieve hashes from kmers/minimizers
            std::vector< size_t > hashes_f;
            std::vector< size_t > hashes_r;

            // hash count
            uint16_t kmers = 0;
            // Best scoring kmer count
            uint16_t max_kmer_count_read = 0;
            if ( read1_len >= wk_size )
            {
                // Count hashes from both pairs if second is given
                if ( read2_len >= wk_size )
                {
                    // FR --> RF <-- reads
                    // Get hashes of the first pair (forward and reverse strand)
                    // get hashed of the seconf pair in the opposite direction
                    get_hashes( rb.seqs[readID],
                                rb.seqs2[readID],
                                hashes_f,
                                hashes_r,
                                hierarchy_config.window_size,
                                hierarchy_config.offset,
                                kmer_hash,
                                minimiser_hash );
                }
                else
                {
                    // single read
                    get_hashes( rb.seqs[readID],
                                hashes_f,
                                hashes_r,
                                hierarchy_config.window_size,
                                hierarchy_config.offset,
                                kmer_hash,
                                minimiser_hash );
                }
                kmers = hashes_f.size();


                // Sum sequence to totals
                if ( hierarchy_first )
                {
                    total.reads_processed++;
                    total.length_processed += read1_len + read2_len;
                }

                // For each filter in the hierarchy
                for ( uint8_t i = 0; i < filters.size(); ++i )
                {


                    // Calculate threshold for cutoff (keep matches above)
                    uint16_t threshold_cutoff = 0;
                    if ( filters[i].filter_config.rel_cutoff >= 0 )
                        threshold_cutoff = threshold_rel( kmers, filters[i].filter_config.rel_cutoff );
                    else if ( filters[i].filter_config.abs_cutoff >= 0 )
                        threshold_cutoff = threshold_abs( kmers,
                                                          filters[i].filter_config.abs_cutoff,
                                                          hierarchy_config.kmer_size,
                                                          hierarchy_config.offset );

                    // reset low threshold_cutoff to just one kmer (0 would match everywhere)
                    if ( threshold_cutoff == 0 )
                        threshold_cutoff = 1;

                    // count and select matches
                    select_matches(
                        filters[i], matches, hashes_f, hashes_r, agents[i], threshold_cutoff, max_kmer_count_read );
                }
            }

            // store read and matches to be printed
            ReadOut read_out( rb.ids[readID] );

            // if read got valid matches (above cutoff)
            if ( max_kmer_count_read > 0 )
            {
                total.reads_classified++;

                // Calculate threshold for filtering (keep matches above)
                uint16_t threshold_filter = get_threshold_filter(
                    hierarchy_config, kmers, max_kmer_count_read, hierarchy_config.kmer_size, hierarchy_config.offset );

                // Filter matches
                uint32_t count_filtered_matches = filter_matches( read_out, matches, rep, threshold_filter );

                if ( run_lca )
                {
                    ReadOut read_out_lca( rb.ids[readID] );
                    if ( count_filtered_matches == 1 )
                    {
                        // just one match, copy read read_out and set as unique
                        read_out_lca = read_out;
                        rep[read_out.matches[0].target].unique_reads++;
                    }
                    else
                    {
                        lca_matches( read_out, read_out_lca, max_kmer_count_read, lca, rep );
                    }

                    if ( config.output_lca )
                        classified_lca_queue.push( read_out_lca );
                }
                else if ( count_filtered_matches == 1 )
                {
                    // Not running lca and has unique match
                    rep[read_out.matches[0].target].unique_reads++;
                }

                if ( config.output_all )
                    classified_all_queue.push( read_out );

                // read classified, continue to the next
                continue;
            }

            // not classified
            if ( !hierarchy_last ) // if there is more levels, store read
            {
                // seqan::appendValue( left_over_reads.ids, rb.ids[readID] );
                // seqan::appendValue( left_over_reads.seqs, rb.seqs[readID] );
                // MOVE?
                left_over_reads.ids.push_back( rb.ids[readID] );
                left_over_reads.seqs.push_back( rb.seqs[readID] );

                if ( rb.paired )
                {
                    // seqan::appendValue( left_over_reads.seqs2, rb.seqs2[readID] );
                    left_over_reads.seqs2.push_back( rb.seqs2[readID] );
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
            pointer_helper->push( left_over_reads );
    }
}

void write_report( TRep& rep, TTax& tax, std::ofstream& out_rep, std::string hierarchy_label, bool run_lca )
{
    for ( auto const& [target, report] : rep )
    {
        if ( report.matches || report.lca_reads || report.unique_reads )
        {
            out_rep << hierarchy_label << '\t' << target << '\t' << report.matches << '\t' << report.unique_reads
                    << '\t' << report.lca_reads;
            if ( run_lca )
            {
                out_rep << '\t' << tax.at( target ).rank << '\t' << tax.at( target ).name;
            }
            out_rep << '\n';
        }
    }
}

size_t load_filter( THIBF& filter, std::string const& input_filter_file )
{
    std::ifstream              is( input_filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );

    uint32_t                                  parsed_version;
    uint64_t                                  window_size;
    seqan3::shape                             shape{};
    uint8_t                                   parts;
    bool                                      compressed;
    std::vector< std::vector< std::string > > bin_path{};

    archive( parsed_version );
    archive( window_size );
    archive( shape );
    archive( parts );
    archive( compressed );
    archive( bin_path );
    archive( filter );

    return filter.user_bins.num_user_bins();
}

size_t load_filter( TIBF& filter, std::string const& input_filter_file )
{
    std::ifstream              is( input_filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( filter );

    return filter.bin_count();
}

TMap load_map( std::string map_file )
{
    TMap          map;
    std::string   line;
    std::ifstream infile;
    infile.open( map_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field );
        // target <tab> binid
        map[std::stoul( fields[1] )] = fields[0];
    }
    infile.close();
    return map;
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
bool load_files( std::vector< Filter< TFilter > >& filters, std::string hierarchy_label, Config& config, bool run_lca )
{
    uint16_t filter_cnt = 0;
    for ( auto const& filter_config : config.parsed_hierarchy[hierarchy_label].filters )
    {
        TMap    map;
        TTax    tax;
        TFilter filter;
        auto    bin_count = load_filter( filter, filter_config.ibf_file );

        if ( !filter_config.map_file.empty() )
        {
            map = load_map( filter_config.map_file );
        }
        else
        {
            // If multiple hierarchies or filters are used, create unique binid target
            std::string binid_prefix = "";
            if ( config.parsed_hierarchy.size() > 1 || config.parsed_hierarchy[hierarchy_label].filters.size() > 1 )
            {
                binid_prefix = hierarchy_label + "-" + std::to_string( filter_cnt ) + "-";
            }
            // fill map with binids as targets
            for ( uint32_t binid = 0; binid < bin_count; ++binid )
            {
                map[binid] = binid_prefix + std::to_string( binid );
            }
        }

        // Check consistency of map file and ibf file
        // map.size can be smaller than bins set on IBF if sequences were removed
        if ( map.size() > bin_count )
        {
            std::cerr << "ERROR: .ibf and .map files are inconsistent." << std::endl;
            return false;
        }

        if ( run_lca )
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
    std::cerr << "ganon-classify    start time: " << timeGanon.begin() << std::endl;
    std::cerr << "loading filters      elapsed: " << timeLoadFilters.elapsed() << " seconds" << std::endl;
    std::cerr << "classifying+printing elapsed: " << timeClassPrint.elapsed() << " seconds" << std::endl;
    std::cerr << "ganon-classify       elapsed: " << timeGanon.elapsed() << " seconds" << std::endl;
    std::cerr << "ganon-classify      end time: " << timeGanon.end() << std::endl;
    std::cerr << std::endl;
}

void print_stats( Stats& stats, const Config& config, const StopClock& timeClassPrint )
{
    const double elapsed_classification = timeClassPrint.elapsed();
    std::cerr << "ganon-classify processed " << stats.total.reads_processed << " sequences ("
              << stats.total.length_processed / 1000000.0 << " Mbp) in " << elapsed_classification << " seconds ("
              << ( stats.total.reads_processed / 1000.0 ) / ( elapsed_classification / 60.0 ) << " Kseq/m, "
              << ( stats.total.length_processed / 1000000.0 ) / ( elapsed_classification / 60.0 ) << " Mbp/m)"
              << std::endl;
    std::cerr << " - " << stats.total.reads_classified << " reads classified ("
              << ( stats.total.reads_classified / static_cast< double >( stats.total.reads_processed ) ) * 100 << "%)"
              << std::endl;
    std::cerr << "   - " << stats.total.unique_matches << " with unique matches ("
              << ( stats.total.unique_matches / static_cast< double >( stats.total.reads_processed ) ) * 100 << "%)"
              << std::endl;
    std::cerr << "   - " << stats.total.reads_classified - stats.total.unique_matches << " with multiple matches ("
              << ( ( stats.total.reads_classified - stats.total.unique_matches )
                   / static_cast< double >( stats.total.reads_processed ) )
                     * 100
              << "%)" << std::endl;

    float avg_matches = stats.total.reads_classified
                            ? ( stats.total.matches / static_cast< double >( stats.total.reads_classified ) )
                            : 0;
    std::cerr << " - " << stats.total.matches << " matches (avg. " << avg_matches << " match/read classified)"
              << std::endl;
    uint64_t total_reads_unclassified = stats.total.reads_processed - stats.total.reads_classified;
    std::cerr << " - " << total_reads_unclassified << " reads unclassified ("
              << ( total_reads_unclassified / static_cast< double >( stats.total.reads_processed ) ) * 100 << "%)"
              << std::endl;

    if ( stats.total.reads_processed < stats.input_reads )
    {
        std::cerr << " - " << stats.input_reads - stats.total.reads_processed
                  << " reads skipped (shorther than k-mer size)" << std::endl;
    }

    if ( config.parsed_hierarchy.size() > 1 )
    {
        std::cerr << std::endl;
        std::cerr << "By database hierarchy:" << std::endl;
        for ( auto const& h : config.parsed_hierarchy )
        {
            std::string hierarchy_label = h.first;
            avg_matches                 = stats.hierarchy_total[hierarchy_label].reads_classified
                                              ? ( stats.hierarchy_total[hierarchy_label].matches
                                  / static_cast< double >( stats.hierarchy_total[hierarchy_label].reads_classified ) )
                                              : 0;
            std::cerr << " - " << hierarchy_label << ": " << stats.hierarchy_total[hierarchy_label].reads_classified
                      << " classified ("
                      << ( stats.hierarchy_total[hierarchy_label].reads_classified
                           / static_cast< double >( stats.total.reads_processed ) )
                             * 100
                      << "%) " << stats.hierarchy_total[hierarchy_label].unique_matches << " unique ("
                      << ( stats.hierarchy_total[hierarchy_label].unique_matches
                           / static_cast< double >( stats.total.reads_processed ) )
                             * 100
                      << "%) "
                      << stats.hierarchy_total[hierarchy_label].reads_classified
                             - stats.hierarchy_total[hierarchy_label].unique_matches
                      << " multiple ("
                      << ( ( stats.hierarchy_total[hierarchy_label].reads_classified
                             - stats.hierarchy_total[hierarchy_label].unique_matches )
                           / static_cast< double >( stats.total.reads_processed ) )
                             * 100
                      << "%) " << stats.hierarchy_total[hierarchy_label].matches << " matches (avg. " << avg_matches
                      << ")" << std::endl;
        }
    }
}

void parse_reads( SafeQueue< ReadBatches >& queue1, Stats& stats, Config const& config )
{
    for ( auto const& reads_file : config.single_reads )
    {
        seqan3::sequence_file_input< raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq > > fin1{
            reads_file
        };
        for ( auto&& rec : fin1 | seqan3::views::chunk( config.n_reads ) )
        {
            ReadBatches rb{ false };
            for ( auto& [id, seq] : rec )
            {
                rb.ids.push_back( std::move( id ) );
                rb.seqs.push_back( std::move( seq ) );
            }
            stats.input_reads += rb.ids.size();
            queue1.push( std::move( rb ) );
        }
    }
    if ( config.paired_reads.size() > 0 )
    {
        for ( uint16_t pair_cnt = 0; pair_cnt < config.paired_reads.size(); pair_cnt += 2 )
        {
            seqan3::sequence_file_input< raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                fin1{ config.paired_reads[pair_cnt] };
            seqan3::sequence_file_input< raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq > >
                fin2{ config.paired_reads[pair_cnt + 1] };
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
                stats.input_reads += rb.ids.size();
                queue1.push( std::move( rb ) );
            }
        }
    }
    queue1.notify_push_over();
}

void write_classified( SafeQueue< ReadOut >& classified_queue, std::ofstream& out )
{
    while ( true )
    {
        ReadOut ro = classified_queue.pop();
        if ( ro.readID != "" )
        {
            for ( uint32_t i = 0; i < ro.matches.size(); ++i )
            {
                out << ro.readID << '\t' << ro.matches[i].target << '\t' << ro.matches[i].kmer_count << '\n';
            }
        }
        else
        {
            break;
        }
    }
}

void write_unclassified( SafeQueue< ReadOut >& unclassified_queue, std::string out_unclassified_file )
{
    std::ofstream out_unclassified( out_unclassified_file );
    while ( true )
    {
        ReadOut rou = unclassified_queue.pop();
        if ( rou.readID != "" )
        {
            out_unclassified << rou.readID << '\n';
        }
        else
        {
            out_unclassified.close();
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
        for ( uint16_t i = 1; i < filters.size(); ++i )
        {
            // merge taxonomies keeping the first one as a default
            merged_tax.insert( filters[i].tax.begin(), filters[i].tax.end() );
        }
        return merged_tax;
    }
}

template < typename TFilter >
void validate_targets_tax( std::vector< Filter< TFilter > > const& filters, TTax& tax, bool quiet )
{
    for ( auto const& filter : filters )
    {
        for ( auto const& [binid, target] : filter.map )
        {
            if ( tax.count( target ) == 0 )
            {
                tax[target] = Node{ "1", "no rank", target };
                if ( !quiet )
                    std::cerr << "WARNING: target [" << target << "] without tax entry, setting parent node to 1 (root)"
                              << std::endl;
            }
        }
    }
}

void pre_process_lca( LCA& lca, TTax& tax )
{
    for ( auto const& [target, node] : tax )
    {
        lca.addEdge( node.parent, target );
    }
    lca.doEulerWalk();
}

} // namespace detail

template < typename TFilter >
bool ganon_classify( Config config )
{

    // Start time count
    StopClock timeGanon;
    timeGanon.start();

    // Initialize variables
    StopClock timeLoadFilters;
    StopClock timeClassPrint;

    detail::Stats stats;
    std::ofstream out_rep; // Set default output stream (file or stdout)
    std::ofstream out_all; // output all file
    std::ofstream out_lca; // output lca file

    // If there's no output prefix, redirect to STDOUT
    if ( config.output_prefix.empty() )
    {
        out_rep.copyfmt( std::cout ); // STDOUT
        out_rep.clear( std::cout.rdstate() );
        out_rep.basic_ios< char >::rdbuf( std::cout.rdbuf() );
    }
    else
    {
        out_rep.open( config.output_prefix + ".rep" );
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

    // Thread for reading input files
    std::future< void > read_task = std::async(
        std::launch::async, detail::parse_reads, std::ref( queue1 ), std::ref( stats ), std::ref( config ) );

    // Thread for printing unclassified reads
    SafeQueue< detail::ReadOut > unclassified_queue;
    std::future< void >          write_unclassified_task;
    if ( config.output_unclassified && !config.output_prefix.empty() )
    {
        write_unclassified_task = std::async( std::launch::async,
                                              detail::write_unclassified,
                                              std::ref( unclassified_queue ),
                                              config.output_prefix + ".unc" );
    }


    // Classify reads iteractively for each hierarchy level
    uint16_t hierarchy_id   = 0;
    uint16_t hierarchy_size = config.parsed_hierarchy.size();
    // Check if tax files are present to run lca
    bool run_lca = config.tax.size();
    for ( auto const& [hierarchy_label, hierarchy_config] : config.parsed_hierarchy )
    {
        ++hierarchy_id;
        bool                                     hierarchy_first = ( hierarchy_id == 1 );
        bool                                     hierarchy_last  = ( hierarchy_id == hierarchy_size );
        std::vector< detail::Filter< TFilter > > filters;
        detail::TTax                             tax;
        LCA                                      lca;

        timeLoadFilters.start();
        bool loaded = detail::load_files( filters, hierarchy_label, config, run_lca );
        if ( !loaded )
        {
            std::cerr << "ERROR: loading ibf, map or tax files" << std::endl;
            return false;
        }
        timeLoadFilters.stop();

        if ( run_lca )
        {
            // merge repeated elements from multiple filters
            tax = detail::merge_tax( filters );
            // if target not found in tax, add node target with parent = "1" (root)
            detail::validate_targets_tax( filters, tax, config.quiet );
            // pre-processing of nodes to generate LCA
            detail::pre_process_lca( lca, tax );
        }


        // Thread for printing classified reads (.lca, .all)
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

        if ( !config.output_prefix.empty() )
        {
            if ( config.output_lca && run_lca )
            {
                if ( hierarchy_first || !config.output_single )
                    out_lca.open( hierarchy_config.output_file_lca );
                else // append if not first and output_single
                    out_lca.open( hierarchy_config.output_file_lca, std::ofstream::app );

                // Start writing thread for lca matches
                write_tasks.emplace_back( std::async( std::launch::async,
                                                      detail::write_classified,
                                                      std::ref( classified_lca_queue ),
                                                      std::ref( out_lca ) ) );
            }
            if ( config.output_all )
            {
                if ( hierarchy_first || !config.output_single )
                    out_all.open( hierarchy_config.output_file_all );
                else // append if not first and output_single
                    out_all.open( hierarchy_config.output_file_all, std::ofstream::app );

                // Start writing thread for all matches
                write_tasks.emplace_back( std::async( std::launch::async,
                                                      detail::write_classified,
                                                      std::ref( classified_all_queue ),
                                                      std::ref( out_all ) ) );
            }
        }

        // One report and total counters for each thread
        std::vector< detail::TRep >  reports( config.threads_classify );
        std::vector< detail::Total > totals( config.threads_classify );

        std::vector< std::future< void > > tasks;
        // Threads for classification
        timeClassPrint.start();
        for ( uint16_t taskNo = 0; taskNo < config.threads_classify; ++taskNo )
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
                                            hierarchy_last,
                                            run_lca ) );
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

        // Sum totals for each thread and report into stats
        stats.add_totals( hierarchy_label, totals );
        stats.add_reports( hierarchy_label, rep );

        // write reports
        detail::write_report( rep, tax, out_rep, hierarchy_label, run_lca );

        // Wait here until all files are written
        for ( auto&& task : write_tasks )
        {
            task.get();
        }
        timeClassPrint.stop();

        // Close file for writing (if not STDOUT)
        if ( !config.output_prefix.empty() )
        {
            if ( config.output_lca )
                out_lca.close();
            if ( config.output_all )
                out_all.close();
        }

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

    out_rep << "#total_classified\t" << stats.total.reads_classified << '\n';
    // account for unclassified and skipped sequences
    out_rep << "#total_unclassified\t" << stats.input_reads - stats.total.reads_classified << '\n';
    if ( !config.output_prefix.empty() )
    {
        out_rep.close();
    }

    timeGanon.stop();

    if ( !config.quiet )
    {
        std::cerr << std::endl;
        if ( config.verbose )
        {
            detail::print_time( timeGanon, timeLoadFilters, timeClassPrint );
        }
        detail::print_stats( stats, config, timeClassPrint );
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
