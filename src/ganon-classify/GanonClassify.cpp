#include "GanonClassify.hpp"

#include <utils/LCA.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <cereal/archives/binary.hpp>


#include <atomic>
#include <cinttypes>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace GanonClassify
{

namespace detail
{

typedef seqan3::interleaved_bloom_filter<> TIbf;

typedef seqan3::interleaved_bloom_filter<>::counting_agent_type< uint16_t > Tagent;

typedef std::unordered_map< std::string, uint16_t > TMatches;

struct Node
{
    std::string parent;
    std::string rank;
    std::string name;
};

struct Rep
{
    std::atomic< uint64_t > direct_matches = 0;
    std::atomic< uint64_t > lca_reads      = 0;
    std::atomic< uint64_t > unique_reads   = 0;
};

typedef std::unordered_map< std::string, Rep >  TRep;
typedef std::unordered_map< std::string, Node > TTax;
typedef std::map< uint32_t, std::string >       TMap;

typedef std::unordered_set< std::string > TUniqueTarget;

struct ReadBatches
{

    ReadBatches()
    {
    }

    ReadBatches( bool _paired )
    {
        paired = _paired;
    }

    ReadBatches( bool _paired, std::vector< std::string > _ids, std::vector< std::vector< seqan3::dna5 > > _seqs )
    {
        paired = _paired;
        ids    = _ids;
        seqs   = _seqs;
    }

    ReadBatches( bool                                       _paired,
                 std::vector< std::string >                 _ids,
                 std::vector< std::vector< seqan3::dna5 > > _seqs,
                 std::vector< std::vector< seqan3::dna5 > > _seqs2 )
    {
        paired = _paired;
        ids    = _ids;
        seqs   = _seqs;
        seqs2  = _seqs2;
    }

    bool                                       paired = false;
    std::vector< std::string >                 ids;
    std::vector< std::vector< seqan3::dna5 > > seqs;
    std::vector< std::vector< seqan3::dna5 > > seqs2;
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

struct Stats
{
    Stats()
    : sumread_len{ 0 }
    , totalReads{ 0 }
    , classifiedReads{ 0 }
    , matches{ 0 }
    {
    }

    std::atomic< uint64_t >           sumread_len;
    uint64_t                          totalReads;
    uint64_t                          classifiedReads;
    uint64_t                          matches;
    std::map< std::string, uint64_t > classifiedReads_hierarchy;
    std::map< std::string, uint64_t > matches_hierarchy;

    void load_rep( std::string hierarchy_label, TRep& rep )
    {
        for ( auto& r : rep )
        {
            classifiedReads_hierarchy[hierarchy_label] += r.second.unique_reads + r.second.lca_reads;
            matches_hierarchy[hierarchy_label] += r.second.direct_matches;
            classifiedReads += r.second.unique_reads + r.second.lca_reads;
            matches += r.second.direct_matches;
        }
    }
};


struct Filter
{
    TIbf         ibf;
    TMap         map;
    TTax         tax;
    FilterConfig filter_config;
};

inline uint16_t get_error( uint16_t read_len, uint8_t kmer_size, uint16_t kmer_count, uint8_t offset )
{
    // Return the optimal number of errors for a certain sequence based on the kmer_count

    // If offset > 1, check weather an extra position is necessary to cover the whole read (last k-mer)
    bool extra_pos = ( read_len - kmer_size ) % offset;

    // - ((offset-1)/kmer_size) adjusts the value to be always below the threshold possible, considering the floor usage
    // to calculate the kmer_count
    return std::ceil( ( read_len - kmer_size + offset * ( -kmer_count + extra_pos + 1 ) )
                          / static_cast< float >( kmer_size )
                      - ( ( offset - 1 ) / static_cast< float >( kmer_size ) ) );
}

inline uint16_t get_threshold_errors( uint16_t read_len, uint8_t kmer_size, uint16_t max_error, uint8_t offset )
{
    // Return threshold (number of kmers) based on an optimal number of errors
    // 1 instead of 0 - meaning that if a higher number of errors are allowed the threshold here is
    // just one kmer match (0 would match every read everywhere)

    // If offset > 1, check weather an extra position is necessary to cover the whole read (last k-mer)
    bool extra_pos = ( read_len - kmer_size ) % offset;

    return read_len + 1u > kmer_size * ( 1u + max_error )
               ? std::floor( ( ( read_len - kmer_size - max_error * kmer_size ) / offset ) + 1 + extra_pos )
               : 1u;
}

inline uint16_t get_threshold_kmers( uint16_t read_len, uint8_t kmer_size, float min_kmers, uint8_t offset )
{
    // Return threshold (number of kmers) based on an percentage of kmers. 0 for anything with at least 1 k-mer
    // ceil -> round-up min # k-mers, floor -> round-down for offset

    // If offset > 1, check weather an extra position is necessary to cover the whole read (last k-mer)
    bool extra_pos = ( read_len - kmer_size ) % offset;

    return min_kmers > 0 ? std::floor( std::ceil( ( read_len - kmer_size ) * min_kmers ) / offset + 1 + extra_pos )
                         : 1u;
}


inline void check_unique( ReadOut& read_out_lca,
                          uint16_t read_len,
                          uint16_t kmer_size,
                          uint16_t max_error_unique,
                          uint8_t  offset,
                          TTax&    tax,
                          TRep&    rep )
{
    uint16_t threshold_error_unique = get_threshold_errors( read_len, kmer_size, max_error_unique, offset );
    // if kmer count is lower than expected set match to parent node
    if ( read_out_lca.matches[0].kmer_count < threshold_error_unique )
    {
        read_out_lca.matches[0].target = tax.at( read_out_lca.matches[0].target ).parent; // parent node
        rep.at( read_out_lca.matches[0].target ).lca_reads++;                             // count as lca for parent
    }
    else
    {
        rep.at( read_out_lca.matches[0].target ).unique_reads++; // count as unique for target
    }
}

void select_matches( TMatches&                matches,
                     std::vector< uint16_t >& selectedBins,
                     std::vector< uint16_t >& selectedBinsRev,
                     Filter&                  filter,
                     uint16_t                 threshold,
                     uint16_t&                maxKmerCountRead )
{
    // for each bin
    // for ( uint32_t binNo = 0; binNo < filter.ibf.noOfBins; ++binNo )
    // loop in map structure to avoid extra validations when map.size() < filter.ibf.noOfBins when ibf is updated and
    // sequences removed also avoid the error of spurius results from empty bins (bug reported)
    for ( auto const& [binNo, target] : filter.map )
    {
        // if kmer count is higher than threshold
        if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
        {
            // get best matching strand
            uint16_t maxKmerCountBin = std::max( selectedBins[binNo], selectedBinsRev[binNo] );
            // keep only the best match target/read when same targets are split in several
            // bins
            if ( matches.count( target ) == 0 || maxKmerCountBin > matches[target] )
            {
                // store match to target
                matches[target] = maxKmerCountBin;
                if ( maxKmerCountBin > maxKmerCountRead )
                    maxKmerCountRead = maxKmerCountBin;
            }
        }
    }
}

uint16_t find_matches( TMatches&                    matches,
                       std::vector< Filter >&       filters,
                       std::vector< Tagent >&       agents,
                       std::vector< seqan3::dna5 >& read_seq,
                       uint8_t                      offset,
                       auto&                        hash_adaptor)
{

    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;

    for(uint8_t i = 0; i < filters.size(); ++i)
    {
        seqan3::counting_vector< uint16_t > selectedBins = agents[i].bulk_count( read_seq | hash_adaptor );
        seqan3::counting_vector< uint16_t > selectedBinsRev =
            agents[i].bulk_count( read_seq | std::views::reverse | seqan3::views::complement | hash_adaptor );

        uint16_t threshold =
            ( filters[i].filter_config.min_kmers > -1 )
                ? get_threshold_kmers(
                      read_seq.size(), filters[i].filter_config.kmer_size, filters[i].filter_config.min_kmers, offset )
                : get_threshold_errors(
                      read_seq.size(), filters[i].filter_config.kmer_size, filters[i].filter_config.max_error, offset );
        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filters[i], threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint16_t find_matches_paired( TMatches&                    matches,
                              std::vector< Filter >&       filters,
                              std::vector< Tagent >&       agents,
                              std::vector< seqan3::dna5 >& read_seq,
                              std::vector< seqan3::dna5 >& read_seq2,
                              uint16_t                     effective_read_len,
                              uint8_t                      offset,
                              auto&                        hash_adaptor )
{
    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;

    for(uint8_t i = 0; i < filters.size(); ++i)
    {

        // IBF count
        // FR
        seqan3::counting_vector< uint16_t > selectedBins = agents[i].bulk_count( read_seq | hash_adaptor );
        selectedBins += agents[i].bulk_count( read_seq2 | std::views::reverse | seqan3::views::complement | hash_adaptor );

        // RF
        seqan3::counting_vector< uint16_t > selectedBinsRev =
            agents[i].bulk_count( read_seq | std::views::reverse | seqan3::views::complement | hash_adaptor );
        selectedBinsRev += agents[i].bulk_count( read_seq2 | hash_adaptor );

        uint16_t threshold =
            ( filters[i].filter_config.min_kmers > -1 )
                ? get_threshold_kmers(
                      effective_read_len, filters[i].filter_config.kmer_size, filters[i].filter_config.min_kmers, offset )
                : get_threshold_errors(
                      effective_read_len, filters[i].filter_config.kmer_size, filters[i].filter_config.max_error, offset );

        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filters[i], threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint32_t filter_matches( ReadOut&  read_out,
                         TMatches& matches,
                         TRep&     rep,
                         uint16_t  len,
                         uint16_t  max_kmer_count_read,
                         uint16_t  kmer_size,
                         uint8_t   offset,
                         int16_t   strata_filter )
{

    uint16_t threshold_strata = 1; // minimum threshold (when strata_filter == -1)
    if ( strata_filter > -1 )
    {
        // get maximum possible number of error for this read
        uint16_t max_error = get_error( len, kmer_size, max_kmer_count_read, offset );
        // get min kmer count necesary to achieve the calculated number of errors
        threshold_strata = get_threshold_errors( len, kmer_size, max_error + strata_filter, offset );
    }

    for ( auto const& v : matches )
    { // matches[target] = kmerCount
        if ( v.second >= threshold_strata )
        { // apply strata filter
            rep.at( v.first ).direct_matches++;
            read_out.matches.push_back( ReadMatch{ v.first, v.second } );
        }
    }

    return read_out.matches.size();
}

void lca_matches( ReadOut& read_out, ReadOut& read_out_lca, uint16_t max_kmer_count_read, LCA& lca, TRep& rep )
{
    std::vector< std::string > targets;
    for ( auto& r : read_out.matches )
    {
        targets.push_back( r.target );
    }

    std::string target_lca = lca.getLCA( targets );
    rep.at( target_lca ).lca_reads++;
    read_out_lca.matches.push_back( ReadMatch{ target_lca, max_kmer_count_read } );
}


void classify( std::vector< Filter >&    filters,
               LCA&                      lca,
               TTax&                     tax,
               TRep&                     rep,
               SafeQueue< ReadOut >&     classified_all_queue,
               SafeQueue< ReadOut >&     classified_lca_queue,
               SafeQueue< ReadOut >&     unclassified_queue,
               Stats&                    stats,
               Config const&             config,
               SafeQueue< ReadBatches >* pointer_current,
               SafeQueue< ReadBatches >* pointer_helper,
               bool                      hierarchy_first,
               bool                      hierarchy_last,
               uint8_t                   offset,
               int16_t                   max_error_unique,
               int16_t                   strata_filter )
{

    // k-mer sizes should be the same among filters
    uint8_t kmer_size = filters[0].filter_config.kmer_size;

    // oner hash adaptor per thread
    auto hash_adaptor = seqan3::views::kmer_hash( seqan3::ungapped{ filters[0].filter_config.kmer_size } );

    // one agent per thread per filter
    std::vector< Tagent > agents;
    for ( Filter& filter : filters )
    {
        agents.push_back( filter.ibf.counting_agent< uint16_t >() );
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
            // receives len of first read
            uint16_t effective_read_len = rb.seqs[readID].size();

            // count lens just once
            if ( hierarchy_first )
                stats.sumread_len += effective_read_len;

            TMatches matches;
            ReadOut  read_out( rb.ids[readID] );
            uint16_t max_kmer_count_read = 0;

            // just skip classification, add read to left over (dbs can have different kmer sizes) or unclassified
            if ( effective_read_len >= kmer_size )
            {
                if ( rb.paired ) // paired-end mode
                {
                    uint16_t read2_len = rb.seqs2[readID].size();
                    if ( read2_len >= kmer_size )
                    {
                        // add to effective length the pair for error calculation
                        effective_read_len += read2_len + 1 - kmer_size;
                        // count lens just once
                        if ( hierarchy_first )
                            stats.sumread_len += read2_len;

                        max_kmer_count_read = find_matches_paired(
                            matches, filters, agents, rb.seqs[readID], rb.seqs2[readID], effective_read_len, offset, hash_adaptor );
                    }
                }
                else // single-end mode
                {
                    max_kmer_count_read = find_matches( matches, filters, agents, rb.seqs[readID], offset, hash_adaptor );
                }
            }

            // if there were some assignments
            if ( max_kmer_count_read > 0 )
            {

                // filter matches
                uint32_t count_filtered_matches = filter_matches( read_out,
                                                                  matches,
                                                                  rep,
                                                                  effective_read_len,
                                                                  max_kmer_count_read,
                                                                  kmer_size,
                                                                  config.offset,
                                                                  strata_filter );

                // If there are matches
                if ( count_filtered_matches > 0 )
                {

                    ReadOut read_out_lca( rb.ids[readID] );
                    if ( count_filtered_matches == 1 )
                    {
                        read_out_lca = read_out; // just one match, copy read read_out
                        if ( max_error_unique >= 0 )
                        {
                            // re-classify read to parent if threshold<=max_error_unique
                            check_unique(
                                read_out_lca, effective_read_len, kmer_size, max_error_unique, config.offset, tax, rep );
                        }
                        else
                        {
                            rep.at( read_out.matches[0].target ).unique_reads++;
                        }
                    }
                    else
                    {
                        lca_matches( read_out, read_out_lca, max_kmer_count_read, lca, rep );
                    }

                    classified_lca_queue.push( read_out_lca );

                    if ( config.output_all )
                        classified_all_queue.push( read_out );

                    // read classified, continue to the next
                    continue;
                }
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

void write_report( TRep&       rep,
                   TTax&       tax,
                   Stats&      stats,
                   std::string output_file_rep,
                   std::string hierarchy_label,
                   bool        hierarchy_first,
                   bool        hierarchy_last )
{
    std::ofstream out_rep;
    if ( hierarchy_first )
        out_rep.open( output_file_rep );
    else // append if not first and output_single
        out_rep.open( output_file_rep, std::ofstream::app );

    for ( auto const& [target, report] : rep )
    {
        if ( report.direct_matches || report.lca_reads || report.unique_reads )
        {
            out_rep << hierarchy_label << '\t' << target << '\t' << report.direct_matches << '\t' << report.unique_reads
                    << '\t' << report.lca_reads << '\t' << tax.at( target ).rank << '\t' << tax.at( target ).name
                    << '\n';
        }
    }

    if ( hierarchy_last )
    {
        out_rep << "#total_classified\t" << stats.classifiedReads << '\n';
        out_rep << "#total_unclassified\t" << stats.totalReads - stats.classifiedReads << '\n';
    }
    out_rep.close();
}

bool load_filters( std::vector< Filter >& filters,
                   TUniqueTarget&         unique_targets,
                   std::string            hierarchy_label,
                   Config&                config )
{
    uint16_t first_kmer_size = 0;
    for ( auto const& filter_config : config.parsed_hierarchy[hierarchy_label].filters )
    {

        TIbf filter;
        TMap map;
        TTax tax;

        // ibf file
        std::ifstream              is( filter_config.ibf_file, std::ios::binary );
        cereal::BinaryInputArchive archive( is );
        archive( filter );

        if ( first_kmer_size == 0 )
        {
            first_kmer_size = filter_config.kmer_size;
        }
        else if ( first_kmer_size != filter_config.kmer_size )
        {
            std::cerr
                << "ERROR: filters on the same hierarchy should have same k-mer size configuration. Ignoring filter: "
                << filter_config.ibf_file << std::endl;
            continue;
        }

        if ( config.offset > first_kmer_size )
        {
            std::cerr << "WARNING: offset cannot be bigger than k-mer size. Setting offset to " << first_kmer_size
                      << std::endl;
            config.offset = first_kmer_size;
        }

        std::string   line;
        std::ifstream infile;

        // map file
        infile.open( filter_config.map_file );
        while ( std::getline( infile, line, '\n' ) )
        {
            std::istringstream         stream_line( line );
            std::vector< std::string > fields;
            std::string                field;
            while ( std::getline( stream_line, field, '\t' ) )
                fields.push_back( field );
            // target <tab> binid
            map[std::stoul( fields[1] )] = fields[0];
            unique_targets.insert( fields[0] );
        }
        infile.close();


        // Check consistency of map file and ibf file
        // map.size can be smaller than bins set on IBF if sequences were removed
        if ( map.size() > filter.bin_count() )
        {
            std::cerr << "ERROR: .ibf and .map files are inconsistent." << std::endl;
            return false;
        }

        // tax file
        infile.open( filter_config.tax_file );
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

        filters.push_back( Filter{ std::move( filter ), map, tax, filter_config } );
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
    std::cerr << "ganon-classify processed " << stats.totalReads << " sequences (" << stats.sumread_len / 1000000.0
              << " Mbp) in " << elapsed_classification << " seconds ("
              << ( stats.totalReads / 1000.0 ) / ( elapsed_classification / 60.0 ) << " Kseq/m, "
              << ( stats.sumread_len / 1000000.0 ) / ( elapsed_classification / 60.0 ) << " Mbp/m)" << std::endl;
    std::cerr << " - " << stats.classifiedReads << " sequences classified ("
              << ( stats.classifiedReads / static_cast< double >( stats.totalReads ) ) * 100 << "%)" << std::endl;

    if ( stats.matches > 0 )
    {
        std::cerr << " - " << stats.matches << " matches (avg. "
                  << ( stats.matches / static_cast< double >( stats.classifiedReads ) ) << " match/read)" << std::endl;
        if ( config.parsed_hierarchy.size() > 1 )
        {
            for ( auto& h : config.parsed_hierarchy )
            {
                std::cerr << "    " << h.first << ": " << stats.classifiedReads_hierarchy[h.first] << " sequences ("
                          << ( stats.classifiedReads_hierarchy[h.first] / static_cast< double >( stats.totalReads ) )
                                 * 100
                          << "%) " << stats.matches_hierarchy[h.first] << " matches (avg. "
                          << ( stats.matches_hierarchy[h.first]
                               / static_cast< double >( stats.classifiedReads_hierarchy[h.first] ) )
                          << ")" << std::endl;
            }
        }
    }

    std::cerr << " - " << stats.totalReads - stats.classifiedReads << " sequences unclassified ("
              << ( ( stats.totalReads - stats.classifiedReads ) / static_cast< double >( stats.totalReads ) ) * 100
              << "%)" << std::endl;
}

/*uint32_t parse_single_reads( seqan::SeqFileIn&                 seqFileIn,
                             uint32_t                          pos_start,
                             uint32_t                          n_reads,
                             SafeQueue< detail::ReadBatches >& queue1 )
{

    seqan::setPosition( seqFileIn, pos_start ); // rewind the file to the last valid position without errors
    seqan::CharString id;
    seqan::CharString seq;

    seqan::StringSet< seqan::CharString > ids;
    seqan::StringSet< seqan::CharString > seqs;
    uint32_t                              read_cnt    = 0;
    uint32_t                              total_reads = 0;
    while ( !seqan::atEnd( seqFileIn ) )
    {
        try
        {
            seqan::readRecord( id, seq, seqFileIn );
            read_cnt++;

            seqan::appendValue( seqs, seq );
            seqan::appendValue( ids, id );
            if ( read_cnt == n_reads )
            {
                queue1.push( detail::ReadBatches{ false, ids, seqs } );
                seqan::clear( ids );
                seqan::clear( seqs );
                total_reads += read_cnt;
                read_cnt = 0;
            }
        }
        catch ( seqan::Exception const& e )
        {
            std::cerr << "ERROR: " << e.what() << " [@" << id << "]" << std::endl;
        }
    }
    // left overs
    if ( seqan::length( ids ) > 0 )
    {
        queue1.push( detail::ReadBatches{ false, ids, seqs } );
        total_reads += read_cnt;
    }

    return total_reads;
}*/

void parse_reads( SafeQueue< detail::ReadBatches >& queue1, Stats& stats, Config const& config )
{
    for ( auto const& reads_file : config.single_reads )
    {
        seqan3::sequence_file_input fin1{ reads_file };
        for ( auto&& rec : fin1 | ranges::views::chunk( config.n_reads ) )
        {
            std::vector< std::string >                 ids;
            std::vector< std::vector< seqan3::dna5 > > seqs;
            for ( auto& [seq, id, qual] : rec )
            {
                ids.push_back( std::move( id ) );
                seqs.push_back( std::move( seq ) );
            }
            stats.totalReads += ids.size();
            queue1.push( detail::ReadBatches{ false, ids, seqs } );
        }
    }
    if ( config.paired_reads.size() > 0 )
    {
        for ( uint16_t pair_cnt = 0; pair_cnt < config.paired_reads.size(); pair_cnt += 2 )
        {

            seqan3::sequence_file_input fin1{ config.paired_reads[pair_cnt] };
            seqan3::sequence_file_input fin2{ config.paired_reads[pair_cnt + 1] };
            for ( auto&& rec : fin1 | ranges::views::chunk( config.n_reads ) )
            {

                // TODO take==chunk?
                auto&& rec2 = fin2 | std::views::take( config.n_reads );

                std::vector< std::string >                 ids;
                std::vector< std::vector< seqan3::dna5 > > seqs;
                std::vector< std::vector< seqan3::dna5 > > seqs2;
                for ( auto& [seq, id, qual] : rec )
                {
                    ids.push_back( std::move( id ) );
                    seqs.push_back( std::move( seq ) );
                }
                for ( auto& [seq, id, qual] : rec2 )
                {
                    seqs2.push_back( std::move( seq ) );
                }
                stats.totalReads += ids.size();
                queue1.push( detail::ReadBatches{ true, ids, seqs, seqs2 } );
            }
        }
    }
    queue1.notify_push_over();
}

void write_classified( SafeQueue< detail::ReadOut >& classified_all_queue, std::ofstream& out )
{
    while ( true )
    {
        detail::ReadOut ro = classified_all_queue.pop();
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

void write_unclassified( SafeQueue< detail::ReadOut >& unclassified_queue, std::string out_unclassified_file )
{
    std::ofstream out_unclassified( out_unclassified_file );
    while ( true )
    {
        detail::ReadOut rou = unclassified_queue.pop();
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

TTax merge_tax( std::vector< detail::Filter >& filters )
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

void validate_targets_tax( const TUniqueTarget unique_targets, TTax& tax )
{
    for ( const auto& target : unique_targets )
    {
        if ( tax.count( target ) == 0 )
        {
            tax[target] = Node{ "1", "no rank", target };
            std::cerr << "WARNING: target [" << target << "] without tax entry, setting parent node to 1 (root)"
                      << std::endl;
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

bool run( Config config )
{

    if ( !config.validate() )
        return false;

    // Time control
    StopClock timeGanon;
    timeGanon.start();
    StopClock timeLoadFilters;
    StopClock timeClassPrint;

    if ( !config.quiet && config.verbose )
    {
        std::cerr << config;
    }

    //////////////////////////////

    // Set output stream (file or stdout)
    std::ofstream out;
    std::ofstream out_all;

    // If there's no output prefix, redirect to STDOUT
    if ( config.output_prefix.empty() )
    {
        out.copyfmt( std::cout ); // STDOUT
        out.clear( std::cout.rdstate() );
        out.basic_ios< char >::rdbuf( std::cout.rdbuf() );
    }

    // Queues for internal read handling
    // queue1 get reads from file
    // queue2 will get unclassified reads if hierachy == 2
    // if hierachy == 3 queue1 is used for unclassified and so on
    SafeQueue< detail::ReadBatches > queue1(
        config.n_batches ); // config.n_batches*config.n_reads = max. amount of reads in memory
    SafeQueue< detail::ReadBatches > queue2;

    // Statistics values
    detail::Stats stats;

    // Thread for reading input files
    std::future< void > read_task = std::async(
        std::launch::async, detail::parse_reads, std::ref( queue1 ), std::ref( stats ), std::ref( config ) );

    // SafeQueue for printing unclassified
    SafeQueue< detail::ReadOut > unclassified_queue;
    // Thread for printing unclassified reads
    std::future< void > write_unclassified_task;
    if ( config.output_unclassified && !config.output_prefix.empty() )
    {
        write_unclassified_task = std::async( std::launch::async,
                                              detail::write_unclassified,
                                              std::ref( unclassified_queue ),
                                              config.output_prefix + ".unc" );
    }


    // Pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_current = &queue1; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_helper  = &queue2; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_extra;             // pointer to the queues

    uint16_t hierarchy_id   = 0;
    uint16_t hierarchy_size = config.parsed_hierarchy.size();

    // For every hierarchy level
    for ( auto const& [hierarchy_label, hierarchy_config] : config.parsed_hierarchy )
    {
        ++hierarchy_id;
        bool hierarchy_first = ( hierarchy_id == 1 );
        bool hierarchy_last  = ( hierarchy_id == hierarchy_size );

        timeLoadFilters.start();
        std::vector< detail::Filter > filters;
        detail::TUniqueTarget         unique_targets;
        bool                          loaded = detail::load_filters( filters, unique_targets, hierarchy_label, config );
        if ( !loaded )
        {
            return false;
        }
        timeLoadFilters.stop();

        // merge repeated elements
        detail::TTax tax = detail::merge_tax( filters );

        // Reports
        detail::TRep rep;

        // initialize entries for the report with all possible taxa (to direct access it in multiple threads)
        for ( const auto& it : tax )
        {
            rep[it.first]; // initialize using [] due to std::atomic being not movable/copyable
        }

        // validate targets and tax
        // if target not found in tax, add node target with parent = "1" (root)
        detail::validate_targets_tax( unique_targets, tax );

        // pre-processing of nodes
        LCA lca;
        detail::pre_process_lca( lca, tax );

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
            if ( hierarchy_first || !config.output_single )
                out.open( hierarchy_config.output_file_lca );
            else // append if not first and output_single
                out.open( hierarchy_config.output_file_lca, std::ofstream::app );

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

        // Start writing thread for lca matches
        write_tasks.emplace_back( std::async(
            std::launch::async, detail::write_classified, std::ref( classified_lca_queue ), std::ref( out ) ) );

        std::vector< std::future< void > > tasks;
        // Threads for classification
        timeClassPrint.start();
        for ( uint16_t taskNo = 0; taskNo < config.threads_classify; ++taskNo )
        {
            tasks.emplace_back( std::async( std::launch::async,
                                            detail::classify,
                                            std::ref( filters ),
                                            std::ref( lca ),
                                            std::ref( tax ),
                                            std::ref( rep ),
                                            std::ref( classified_all_queue ),
                                            std::ref( classified_lca_queue ),
                                            std::ref( unclassified_queue ),
                                            std::ref( stats ),
                                            std::ref( config ),
                                            pointer_current,
                                            pointer_helper,
                                            hierarchy_first,
                                            hierarchy_last,
                                            config.offset,
                                            hierarchy_config.max_error_unique,
                                            hierarchy_config.strata_filter ) );
        }

        // Wait here until classification is over
        for ( auto&& task : tasks )
        {
            task.get();
        }

        // Inform that no more reads are going to be pushed
        classified_all_queue.notify_push_over();
        classified_lca_queue.notify_push_over();

        // load rep data into stats
        stats.load_rep( hierarchy_label, rep );

        // write reports
        if ( !config.output_prefix.empty() )
        {
            detail::write_report(
                rep, tax, stats, hierarchy_config.output_file_rep, hierarchy_label, hierarchy_first, hierarchy_last );
        }
        // Wait here until all files are written
        for ( auto&& task : write_tasks )
        {
            task.get();
        }
        timeClassPrint.stop();

        // Close file for writing (if not STDOUT)
        if ( !config.output_prefix.empty() )
        {
            out.close();
            if ( config.output_all )
            {
                out_all.close();
            }
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

} // namespace GanonClassify
