#include "GanonClassify.hpp"

#include <utils/LCA.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <cinttypes>
#include <cmath>
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

typedef seqan3::interleaved_bloom_filter<>                                  TFilter;
typedef seqan3::interleaved_bloom_filter<>::counting_agent_type< uint16_t > TAgent;
typedef std::unordered_map< std::string, uint16_t >                         TMatches;

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

struct Rep
{
    uint64_t matches      = 0;
    uint64_t lca_reads    = 0;
    uint64_t unique_reads = 0;
};

typedef std::unordered_map< std::string, Rep >  TRep;
typedef std::unordered_map< std::string, Node > TTax;
typedef std::map< uint32_t, std::string >       TMap;

struct Total
{
    uint64_t reads_processed  = 0;
    uint64_t length_processed = 0;
    uint64_t reads_classified = 0;
    uint64_t matches          = 0;
    uint64_t unique_reads     = 0;
};

struct Stats
{
    Total total;
    // sum of reports by each hiearchy
    std::map< std::string, Rep > report_summary;
    // reads classified for each hiearchy
    std::map< std::string, uint64_t > reads_classified;
    // unique reads for each hiearchy
    std::map< std::string, uint64_t > unique_reads;

    void add_totals( std::string hierarchy_label, std::vector< Total > const& totals )
    {
        // add several totals (from threads)into one
        for ( auto const& t : totals )
        {
            // total.reads_processed+=t.reads_processed; // filled while reading files
            total.length_processed += t.length_processed;
            total.reads_classified += t.reads_classified;
            reads_classified[hierarchy_label] += t.reads_classified;
        }
    }

    void add_report_summary( std::string hierarchy_label, TRep const& report )
    {
        // sum all report numbers for each target to the hiearchy
        for ( auto const& [target, rep] : report )
        {
            total.matches += rep.matches;
            total.unique_reads += rep.unique_reads;
            unique_reads[hierarchy_label] += rep.unique_reads;
            report_summary[hierarchy_label].matches += rep.matches;
            report_summary[hierarchy_label].lca_reads += rep.lca_reads;
            report_summary[hierarchy_label].unique_reads += rep.unique_reads;
        }
    }
};

struct Filter
{
    TFilter      ibf;
    TMap         map;
    TTax         tax;
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


inline uint16_t get_error(
    uint16_t read1_len, uint16_t read2_len, uint8_t kmer_size, uint16_t kmer_count, uint8_t offset )
{
    // Return the optimal number of errors for a certain sequence based on the kmer_count

    // If offset > 1, check weather an extra position is necessary to cover the whole read (last k-mer)
    bool extra_pos = ( read1_len - kmer_size ) % offset;

    // If second read is present, account all possible kmers to equation
    uint16_t read2_max_kmers = 0;
    if ( read2_len )
    {
        read2_max_kmers = ( ( read2_len - kmer_size ) / offset ) + ( ( read1_len - kmer_size ) % offset ) + 1;
    }

    // - ((offset-1)/kmer_size) adjusts the value to be always below the threshold possible, considering the floor usage
    // to calculate the kmer_count
    return std::ceil( ( read1_len - kmer_size + offset * ( -kmer_count + read2_max_kmers + extra_pos + 1 ) )
                          / static_cast< float >( kmer_size )
                      - ( ( offset - 1 ) / static_cast< float >( kmer_size ) ) );
}

inline int16_t get_threshold_errors( uint16_t read_len, uint8_t kmer_size, uint16_t max_error, uint8_t offset )
{
    // Return threshold (number of kmers) based on an optimal number of errors

    // If offset > 1, check weather an extra position is necessary to cover the whole read (last k-mer)
    bool extra_pos = ( read_len - kmer_size ) % offset;

    return std::floor( ( ( read_len - kmer_size - max_error * kmer_size ) / offset ) + 1 + extra_pos );
}

inline int16_t get_threshold_kmers( uint16_t read_len, uint8_t kmer_size, float min_kmers, uint8_t offset )
{
    // Return threshold (number of kmers) based on an percentage of kmers
    // ceil -> round-up min # k-mers, floor -> round-down for offset

    // If offset > 1, check weather an extra position is necessary to cover the whole read (last k-mer)
    bool extra_pos = ( read_len - kmer_size ) % offset;

    // min_kmers==0 return everything with at least one kmer match
    return min_kmers > 0 ? std::floor( std::ceil( ( read_len - kmer_size ) * min_kmers ) / offset + 1 + extra_pos )
                         : 1u;
}


inline void check_unique( ReadOut& read_out_lca,
                          uint16_t read1_len,
                          uint16_t read2_len,
                          uint8_t  kmer_size,
                          uint16_t max_error_unique,
                          uint8_t  offset,
                          TTax&    tax,
                          TRep&    rep )
{
    int16_t threshold_error_unique = get_threshold_errors( read1_len, kmer_size, max_error_unique, offset );
    if ( read2_len )
    {
        threshold_error_unique += get_threshold_errors( read2_len, kmer_size, max_error_unique, offset );
    }
    // if kmer count is lower than expected set match to parent node
    if ( read_out_lca.matches[0].kmer_count < threshold_error_unique )
    {
        read_out_lca.matches[0].target = tax.at( read_out_lca.matches[0].target ).parent; // parent node
        rep[read_out_lca.matches[0].target].lca_reads++;                                  // count as lca for parent
    }
    else
    {
        rep[read_out_lca.matches[0].target].unique_reads++; // count as unique for target
    }
}

void select_matches( TMatches&                matches,
                     std::vector< uint16_t >& selectedBins,
                     std::vector< uint16_t >& selectedBinsRev,
                     Filter&                  filter,
                     int16_t                  threshold,
                     uint16_t&                maxKmerCountRead )
{
    // Threshold can be negative. If a higher number of errors are allowed the threshold here is
    // just one kmer match (0 would match every read everywhere)
    if ( threshold < 1 )
        threshold = 1;

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

auto get_offset_hashes( auto& hashes, uint8_t offset )
{
    // return offset of hashes, always keep last one to cover full read (extra_pos)
    int  n_kmers     = hashes.size() - 1;
    auto nth_element = [offset, n_kmers]( auto&& tuple ) {
        return ( std::get< 0 >( tuple ) % offset == 0 || n_kmers == std::get< 0 >( tuple ) );
    };
    return seqan3::views::zip( std::views::iota( 0 ), hashes ) | std::views::filter( nth_element ) | std::views::values;
}

uint16_t find_matches( TMatches&                    matches,
                       std::vector< Filter >&       filters,
                       std::vector< TAgent >&       agents,
                       std::vector< seqan3::dna5 >& read_seq,
                       uint8_t                      kmer_size,
                       uint8_t                      offset,
                       auto&                        hash_adaptor )
{


    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;

    for ( uint8_t i = 0; i < filters.size(); ++i )
    {
        seqan3::counting_vector< uint16_t > selectedBins;
        seqan3::counting_vector< uint16_t > selectedBinsRev;

        if ( offset > 1 )
        {
            auto hashes{ read_seq | hash_adaptor | seqan3::views::to< std::vector > };
            auto hashesRev{ read_seq | std::views::reverse | seqan3::views::complement | hash_adaptor
                            | seqan3::views::to< std::vector > };

            selectedBins    = agents[i].bulk_count( get_offset_hashes( hashes, offset ) );
            selectedBinsRev = agents[i].bulk_count( get_offset_hashes( hashesRev, offset ) );
        }
        else
        {
            selectedBins = agents[i].bulk_count( read_seq | hash_adaptor );
            selectedBinsRev =
                agents[i].bulk_count( read_seq | std::views::reverse | seqan3::views::complement | hash_adaptor );
        }

        int16_t threshold =
            ( filters[i].filter_config.min_kmers > -1 )
                ? get_threshold_kmers( read_seq.size(), kmer_size, filters[i].filter_config.min_kmers, offset )
                : get_threshold_errors( read_seq.size(), kmer_size, filters[i].filter_config.max_error, offset );

        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filters[i], threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint16_t find_matches_paired( TMatches&                    matches,
                              std::vector< Filter >&       filters,
                              std::vector< TAgent >&       agents,
                              std::vector< seqan3::dna5 >& read_seq,
                              std::vector< seqan3::dna5 >& read_seq2,
                              uint8_t                      kmer_size,
                              uint8_t                      offset,
                              auto&                        hash_adaptor )
{

    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;

    for ( uint8_t i = 0; i < filters.size(); ++i )
    {

        seqan3::counting_vector< uint16_t > selectedBins;
        seqan3::counting_vector< uint16_t > selectedBinsRev;

        if ( offset > 1 )
        {
            auto hashes{ read_seq | hash_adaptor | seqan3::views::to< std::vector > };
            auto hashes2{ read_seq2 | hash_adaptor | seqan3::views::to< std::vector > };
            auto hashesRev{ read_seq | std::views::reverse | seqan3::views::complement | hash_adaptor
                            | seqan3::views::to< std::vector > };
            auto hashes2Rev{ read_seq2 | std::views::reverse | seqan3::views::complement | hash_adaptor
                             | seqan3::views::to< std::vector > };

            selectedBins = agents[i].bulk_count( get_offset_hashes( hashes, offset ) );
            selectedBins += agents[i].bulk_count( get_offset_hashes( hashes2Rev, offset ) );

            selectedBinsRev = agents[i].bulk_count( get_offset_hashes( hashesRev, offset ) );
            selectedBinsRev += agents[i].bulk_count( get_offset_hashes( hashes2, offset ) );
        }
        else
        {
            // FR
            selectedBins = agents[i].bulk_count( read_seq | hash_adaptor );
            selectedBins +=
                agents[i].bulk_count( read_seq2 | std::views::reverse | seqan3::views::complement | hash_adaptor );

            // RF
            selectedBinsRev =
                agents[i].bulk_count( read_seq | std::views::reverse | seqan3::views::complement | hash_adaptor );
            selectedBinsRev += agents[i].bulk_count( read_seq2 | hash_adaptor );
        }

        // Calculate error rate on one read
        int16_t threshold =
            ( filters[i].filter_config.min_kmers > -1 )
                ? get_threshold_kmers( read_seq.size(), kmer_size, filters[i].filter_config.min_kmers, offset )
                : get_threshold_errors( read_seq.size(), kmer_size, filters[i].filter_config.max_error, offset );

        // sum kmers of second read (if using errors, get all possible kmers)
        threshold +=
            ( filters[i].filter_config.min_kmers > -1 )
                ? get_threshold_kmers( read_seq2.size(), kmer_size, filters[i].filter_config.min_kmers, offset )
                : get_threshold_errors( read_seq2.size(), kmer_size, 0, offset );

        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filters[i], threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint32_t filter_matches( ReadOut&  read_out,
                         TMatches& matches,
                         TRep&     rep,
                         uint16_t  read1_len,
                         uint16_t  read2_len,
                         uint16_t  max_kmer_count_read,
                         uint8_t   kmer_size,
                         uint8_t   offset,
                         int16_t   strata_filter )
{

    int16_t threshold_strata = 1; // minimum threshold (when strata_filter == -1)
    if ( strata_filter > -1 )
    {
        // get maximum possible number of errors
        uint16_t max_error = get_error( read1_len, read2_len, kmer_size, max_kmer_count_read, offset );

        // get min kmer count necesary to achieve the calculated number of errors
        threshold_strata = get_threshold_errors( read1_len, kmer_size, max_error + strata_filter, offset );

        if ( read2_len )
        {
            threshold_strata += get_threshold_errors( read2_len, kmer_size, 0, offset );
        }
    }

    for ( auto const& v : matches )
    { // matches[target] = kmerCount
        if ( v.second >= threshold_strata )
        { // apply strata filter
            rep[v.first].matches++;
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
    rep[target_lca].lca_reads++;
    read_out_lca.matches.push_back( ReadMatch{ target_lca, max_kmer_count_read } );
}


void classify( std::vector< Filter >&    filters,
               LCA&                      lca,
               TTax&                     tax,
               TRep&                     rep,
               Total&                    total,
               SafeQueue< ReadOut >&     classified_all_queue,
               SafeQueue< ReadOut >&     classified_lca_queue,
               SafeQueue< ReadOut >&     unclassified_queue,
               Config const&             config,
               SafeQueue< ReadBatches >* pointer_current,
               SafeQueue< ReadBatches >* pointer_helper,
               bool                      hierarchy_first,
               bool                      hierarchy_last,
               uint8_t                   offset,
               int16_t                   max_error_unique,
               int8_t                    kmer_size,
               int16_t                   strata_filter,
               bool                      run_lca )
{

    // oner hash adaptor per thread
    auto hash_adaptor = seqan3::views::kmer_hash( seqan3::ungapped{ kmer_size } );

    // one agent per thread per filter
    std::vector< TAgent > agents;
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
            uint16_t read1_len = rb.seqs[readID].size();
            uint16_t read2_len = 0;

            // count lens just once
            if ( hierarchy_first )
                total.length_processed += read1_len;

            TMatches matches;
            ReadOut  read_out( rb.ids[readID] );
            uint16_t max_kmer_count_read = 0;

            // just skip classification, add read to left over (dbs can have different kmer sizes) or unclassified
            if ( read1_len >= kmer_size )
            {
                if ( rb.paired ) // paired-end mode
                {
                    read2_len = rb.seqs2[readID].size();
                    if ( read2_len >= kmer_size )
                    {

                        // count lens just once
                        if ( hierarchy_first )
                            total.length_processed += read2_len;

                        max_kmer_count_read = find_matches_paired( matches,
                                                                   filters,
                                                                   agents,
                                                                   rb.seqs[readID],
                                                                   rb.seqs2[readID],
                                                                   kmer_size,
                                                                   offset,
                                                                   hash_adaptor );
                    }
                }
                else // single-end mode
                {
                    max_kmer_count_read =
                        find_matches( matches, filters, agents, rb.seqs[readID], kmer_size, offset, hash_adaptor );
                }
            }

            // if there were some assignments
            if ( max_kmer_count_read > 0 )
            {

                // filter matches
                uint32_t count_filtered_matches = filter_matches( read_out,
                                                                  matches,
                                                                  rep,
                                                                  read1_len,
                                                                  read2_len,
                                                                  max_kmer_count_read,
                                                                  kmer_size,
                                                                  config.offset,
                                                                  strata_filter );

                // If there are valid matches after filtering
                if ( count_filtered_matches > 0 )
                {
                    total.reads_classified++;

                    if ( run_lca )
                    {
                        ReadOut read_out_lca( rb.ids[readID] );
                        if ( count_filtered_matches == 1 )
                        {
                            read_out_lca = read_out; // just one match, copy read read_out
                            if ( max_error_unique >= 0 )
                            {
                                // re-classify read to parent if threshold<=max_error_unique
                                check_unique( read_out_lca,
                                              read1_len,
                                              read2_len,
                                              kmer_size,
                                              max_error_unique,
                                              config.offset,
                                              tax,
                                              rep );
                            }
                            else
                            {
                                rep[read_out.matches[0].target].unique_reads++;
                            }
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
                        // Not running lca and unique match
                        rep[read_out.matches[0].target].unique_reads++;
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

TFilter load_filter( std::string const& input_filter_file )
{

    TFilter                    filter;
    std::ifstream              is( input_filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( filter );
    return filter;
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

bool load_files( std::vector< Filter >& filters, std::string hierarchy_label, Config& config, bool run_lca )
{
    for ( auto const& filter_config : config.parsed_hierarchy[hierarchy_label].filters )
    {

        TMap    map;
        TTax    tax;
        TFilter filter = load_filter( filter_config.ibf_file );

        if ( !filter_config.map_file.empty() )
        {
            map = load_map( filter_config.map_file );
        }
        else
        {
            // fill map with binids as targets
            for ( uint32_t binid = 0; binid < filter.bin_count(); ++binid )
            {
                map[binid] = std::to_string( binid );
            }
        }

        // Check consistency of map file and ibf file
        // map.size can be smaller than bins set on IBF if sequences were removed
        if ( map.size() > filter.bin_count() )
        {
            std::cerr << "ERROR: .ibf and .map files are inconsistent." << std::endl;
            return false;
        }

        if ( run_lca )
            tax = load_tax( filter_config.tax_file );

        filters.push_back( Filter{ std::move( filter ), std::move( map ), std::move( tax ), filter_config } );
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
    std::cerr << "   - " << stats.total.unique_reads << " with unique matches ("
              << ( stats.total.unique_reads / static_cast< double >( stats.total.reads_processed ) ) * 100 << "%)"
              << std::endl;
    std::cerr << "   - " << stats.total.reads_classified - stats.total.unique_reads << " with multiple matches ("
              << ( ( stats.total.reads_classified - stats.total.unique_reads )
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

    if ( config.parsed_hierarchy.size() > 1 )
    {
        std::cerr << std::endl;
        std::cerr << "By database hierarchy:" << std::endl;
        for ( auto const& h : config.parsed_hierarchy )
        {
            std::string hierarchy_label = h.first;
            avg_matches                 = stats.reads_classified[hierarchy_label]
                              ? ( stats.report_summary[hierarchy_label].matches
                                  / static_cast< double >( stats.reads_classified[hierarchy_label] ) )
                              : 0;
            std::cerr << " - " << hierarchy_label << ": " << stats.reads_classified[hierarchy_label] << " classified ("
                      << ( stats.reads_classified[hierarchy_label]
                           / static_cast< double >( stats.total.reads_processed ) )
                             * 100
                      << "%) " << stats.unique_reads[hierarchy_label] << " unique ("
                      << ( stats.unique_reads[hierarchy_label] / static_cast< double >( stats.total.reads_processed ) )
                             * 100
                      << "%) " << stats.reads_classified[hierarchy_label] - stats.unique_reads[hierarchy_label]
                      << " multiple ("
                      << ( ( stats.reads_classified[hierarchy_label] - stats.unique_reads[hierarchy_label] )
                           / static_cast< double >( stats.total.reads_processed ) )
                             * 100
                      << "%) " << stats.report_summary[hierarchy_label].matches << " matches (avg. " << avg_matches
                      << ")" << std::endl;
        }
    }
}

void parse_reads( SafeQueue< detail::ReadBatches >& queue1, Stats& stats, Config const& config )
{
    for ( auto const& reads_file : config.single_reads )
    {
        seqan3::sequence_file_input fin1{ reads_file };
        for ( auto&& rec : fin1 | seqan3::views::chunk( config.n_reads ) )
        {
            detail::ReadBatches rb{ false };
            for ( auto& [seq, id, qual] : rec )
            {
                rb.ids.push_back( std::move( id ) );
                rb.seqs.push_back( std::move( seq ) );
            }
            stats.total.reads_processed += rb.ids.size();
            queue1.push( std::move( rb ) );
        }
    }
    if ( config.paired_reads.size() > 0 )
    {
        for ( uint16_t pair_cnt = 0; pair_cnt < config.paired_reads.size(); pair_cnt += 2 )
        {
            seqan3::sequence_file_input fin1{ config.paired_reads[pair_cnt] };
            seqan3::sequence_file_input fin2{ config.paired_reads[pair_cnt + 1] };
            for ( auto&& rec : fin1 | seqan3::views::chunk( config.n_reads ) )
            {
                detail::ReadBatches rb{ true };
                for ( auto& [seq, id, qual] : rec )
                {
                    rb.ids.push_back( std::move( id ) );
                    rb.seqs.push_back( std::move( seq ) );
                }
                // loop in the second file and get same amount of reads
                for ( auto& [seq, id, qual] : fin2 | std::views::take( config.n_reads ) )
                {
                    rb.seqs2.push_back( std::move( seq ) );
                }
                stats.total.reads_processed += rb.ids.size();
                queue1.push( std::move( rb ) );
            }
        }
    }
    queue1.notify_push_over();
}

void write_classified( SafeQueue< detail::ReadOut >& classified_queue, std::ofstream& out )
{
    while ( true )
    {
        detail::ReadOut ro = classified_queue.pop();
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

TTax merge_tax( std::vector< detail::Filter > const& filters )
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

void validate_targets_tax( std::vector< detail::Filter > const& filters, TTax& tax )
{
    for ( auto const& filter : filters )
    {
        for ( auto const& [binid, target] : filter.map )
        {
            if ( tax.count( target ) == 0 )
            {
                tax[target] = Node{ "1", "no rank", target };
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

bool run( Config config )
{

    // Validate configuration input
    if ( !config.validate() )
        return false;

    // Print config
    if ( !config.quiet && config.verbose )
        std::cerr << config;

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
        bool                          hierarchy_first = ( hierarchy_id == 1 );
        bool                          hierarchy_last  = ( hierarchy_id == hierarchy_size );
        std::vector< detail::Filter > filters;
        detail::TTax                  tax;
        LCA                           lca;

        timeLoadFilters.start();
        bool loaded = detail::load_files( filters, hierarchy_label, config, run_lca );
        if ( !loaded )
            return false;
        timeLoadFilters.stop();

        if ( run_lca )
        {
            // merge repeated elements from multiple filters
            tax = detail::merge_tax( filters );
            // if target not found in tax, add node target with parent = "1" (root)
            detail::validate_targets_tax( filters, tax );
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
                                            detail::classify,
                                            std::ref( filters ),
                                            std::ref( lca ),
                                            std::ref( tax ),
                                            std::ref( reports[taskNo] ),
                                            std::ref( totals[taskNo] ),
                                            std::ref( classified_all_queue ),
                                            std::ref( classified_lca_queue ),
                                            std::ref( unclassified_queue ),
                                            std::ref( config ),
                                            pointer_current,
                                            pointer_helper,
                                            hierarchy_first,
                                            hierarchy_last,
                                            config.offset,
                                            hierarchy_config.max_error_unique,
                                            hierarchy_config.kmer_size,
                                            hierarchy_config.strata_filter,
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
        stats.add_report_summary( hierarchy_label, rep );

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
    out_rep << "#total_unclassified\t" << stats.total.reads_processed - stats.total.reads_classified << '\n';
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

} // namespace GanonClassify
