#include "GanonClassify.hpp"

#include <utils/LCA.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <seqan/binning_directory.h>

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
#include <vector>

namespace GanonClassify
{

namespace detail
{

// Filter is created with seqan::Offset<1>
// when using seqan::count it acts like seqan::Normal if offset=1
#ifdef GANON_OFFSET
typedef seqan::Offset< 1 > THashCount;
#else
typedef seqan::Normal THashCount;
#endif

typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                 seqan::BDConfig< seqan::Dna5, THashCount, seqan::Uncompressed > >
    TIbf;

typedef seqan::ModifiedString< seqan::ModifiedString< seqan::Dna5String, seqan::ModComplementDna >, seqan::ModReverse >
    TSeqRevComp;

typedef std::unordered_map< std::string, int16_t > TMatches;

struct Node
{
    std::string parent;
    std::string rank;
    std::string name;
};

struct Rep
{
    std::atomic< uint64_t > direct_matches;
    std::atomic< uint64_t > lca_matches;
    std::atomic< uint64_t > unique_matches;
};

typedef std::unordered_map< std::string, Rep >  TRep;
typedef std::unordered_map< std::string, Node > TTax;
typedef std::map< uint32_t, std::string >       TMap;

struct ReadBatches
{

    ReadBatches()
    {
    }

    ReadBatches( bool _paired )
    {
        paired = _paired;
    }

    ReadBatches( bool _paired, seqan::StringSet< seqan::CharString > _ids, seqan::StringSet< seqan::Dna5String > _seqs )
    {
        paired = _paired;
        ids    = _ids;
        seqs   = _seqs;
    }

    ReadBatches( bool                                  _paired,
                 seqan::StringSet< seqan::CharString > _ids,
                 seqan::StringSet< seqan::Dna5String > _seqs,
                 seqan::StringSet< seqan::Dna5String > _seqs2 )
    {
        paired = _paired;
        ids    = _ids;
        seqs   = _seqs;
        seqs2  = _seqs2;
    }

    bool                                  paired = false;
    seqan::StringSet< seqan::CharString > ids;
    seqan::StringSet< seqan::Dna5String > seqs;
    seqan::StringSet< seqan::Dna5String > seqs2;
};

struct ReadMatch
{
    std::string target;
    int16_t     kmer_count;
};

struct ReadOut
{
    ReadOut()
    {
    }

    ReadOut( seqan::CharString _readID )
    {
        readID = _readID;
    }

    seqan::CharString        readID;
    std::vector< ReadMatch > matches;
};

struct Stats
{
    Stats()
    : sumReadLen{ 0 }
    , classifiedReads{ 0 }
    , matches{ 0 }
    , totalReads{ 0 }
    {
    }

    std::atomic< uint64_t > sumReadLen;
    std::atomic< uint64_t > classifiedReads;
    std::atomic< uint64_t > matches;
    uint64_t                totalReads;
};


struct Filter
{
    TIbf         ibf;
    TMap         map;
    TTax         tax;
    FilterConfig filter_config;
};

inline uint16_t get_error( uint16_t readLen, uint16_t kmerSize, uint16_t kmer_count, uint16_t offset )
{
    // Return the optimal number of errors for a certain sequence based on the kmer_count
    // (offset-1) -> to correct for the floor left overs
    return std::ceil( ( -kmerSize + readLen - ( kmer_count * offset + ( offset - 1 ) ) + 1 )
                      / static_cast< float >( kmerSize ) );
}

inline uint16_t get_threshold_errors( uint16_t readLen, uint16_t kmerSize, uint16_t max_error, uint16_t offset )
{
    // Return threshold (number of kmers) based on an optimal number of errors
    // 1 instead of 0 - meaning that if a higher number of errors are allowed the threshold here is
    // just one kmer match (0 would match every read everywhere)
    return readLen + 1u > kmerSize * ( 1u + max_error )
               ? std::floor( ( readLen - kmerSize + 1u - ( max_error * kmerSize ) ) / offset )
               : 1u;
}

inline uint16_t get_threshold_kmers( uint16_t readLen, uint16_t kmerSize, float min_kmers, uint16_t offset )
{
    // Return threshold (number of kmers) based on an percentage of kmers. 0 for anything with at least 1 k-mer
    // ceil -> round-up min # k-mers, floor -> round-down for offset
    return min_kmers > 0 ? std::floor( std::ceil( ( readLen - kmerSize + 1 ) * min_kmers ) / offset ) : 1u;
}


inline bool check_unique(
    ReadOut& read_out, uint16_t read_len, uint16_t kmer_size, uint16_t max_error_unique, uint16_t offset, TTax& tax )
{
    uint16_t threshold_error_unique = get_threshold_errors( read_len, kmer_size, max_error_unique, offset );
    // if kmer count is lower than expected set match to parent node
    if ( read_out.matches[0].kmer_count < threshold_error_unique )
    {
        read_out.matches[0].target = tax[read_out.matches[0].target].parent; // parent node
        return false;
    }
    return true;
}

void select_matches( TMatches&                matches,
                     std::vector< uint16_t >& selectedBins,
                     std::vector< uint16_t >& selectedBinsRev,
                     Filter&                  filter,
                     uint16_t                 threshold,
                     uint16_t&                maxKmerCountRead )
{
    // for each bin
    for ( uint32_t binNo = 0; binNo < filter.ibf.noOfBins; ++binNo )
    {
        // if kmer count is higher than threshold
        if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
        {
            // get best matching strand
            uint16_t maxKmerCountBin = std::max( selectedBins[binNo], selectedBinsRev[binNo] );
            // keep only the best match target/read when same targets are split in several
            // bins
            if ( matches.count( filter.map[binNo] ) == 0 || maxKmerCountBin > matches[filter.map[binNo]] )
            {
                // store match to target
                matches[filter.map[binNo]] = maxKmerCountBin;
                if ( maxKmerCountBin > maxKmerCountRead )
                    maxKmerCountRead = maxKmerCountBin;
            }
        }
    }
}

uint16_t find_matches( TMatches& matches, std::vector< Filter >& filters, seqan::Dna5String& read_seq )
{
    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;
    for ( Filter& filter : filters )
    {
        // IBF count
        std::vector< uint16_t > selectedBins    = seqan::count( filter.ibf, read_seq );
        std::vector< uint16_t > selectedBinsRev = seqan::count( filter.ibf, TSeqRevComp( read_seq ) );

        uint16_t threshold = ( filter.filter_config.min_kmers > 0 )
                                 ? get_threshold_kmers( seqan::length( read_seq ),
                                                        filter.ibf.kmerSize,
                                                        filter.filter_config.min_kmers,
                                                        filter.ibf.offset )
                                 : get_threshold_errors( seqan::length( read_seq ),
                                                         filter.ibf.kmerSize,
                                                         filter.filter_config.max_error,
                                                         filter.ibf.offset );

        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filter, threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint16_t find_matches_paired( TMatches&              matches,
                              std::vector< Filter >& filters,
                              seqan::Dna5String&     read_seq,
                              seqan::Dna5String&     read_seq2,
                              uint16_t               effective_read_len )
{
    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;
    for ( Filter& filter : filters )
    {
        // IBF count
        // FR
        std::vector< uint16_t > selectedBins = seqan::count( filter.ibf, read_seq );
        filter.ibf.count< THashCount >( selectedBins, TSeqRevComp( read_seq2 ) );
        // RF
        std::vector< uint16_t > selectedBinsRev = seqan::count( filter.ibf, TSeqRevComp( read_seq ) );
        filter.ibf.count< THashCount >( selectedBinsRev, read_seq2 );

        uint16_t threshold =
            ( filter.filter_config.min_kmers > 0 )
                ? get_threshold_kmers(
                      effective_read_len, filter.ibf.kmerSize, filter.filter_config.min_kmers, filter.ibf.offset )
                : get_threshold_errors(
                      effective_read_len, filter.ibf.kmerSize, filter.filter_config.max_error, filter.ibf.offset );

        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filter, threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint32_t filter_matches( ReadOut&  read_out,
                         TMatches& matches,
                         TRep&     rep,
                         uint16_t  len,
                         uint16_t  max_kmer_count_read,
                         uint16_t  kmer_size,
                         uint16_t  offset )
{

    // get maximum possible number of error for this read
    uint16_t max_error = get_error( len, kmer_size, max_kmer_count_read, offset );
    // get min kmer count necesary to achieve the calculated number of errors
    uint16_t threshold_strata = get_threshold_errors( len, kmer_size, max_error, offset );

    for ( auto const& v : matches )
    { // matches[target] = kmerCount
        if ( v.second >= threshold_strata )
        { // apply strata filter
            rep[v.first].direct_matches++;
            read_out.matches.push_back( ReadMatch{ v.first, v.second } );
        }
    }

    return read_out.matches.size();
}

void lca_matches( ReadOut& read_out, ReadOut& read_out_lca, LCA& lca, TRep& rep )
{
    std::vector< std::string > targets;
    int16_t                    max_kmer_count = 0;
    for ( auto const & [ target, kmer_count ] : read_out.matches )
    {
        targets.push_back( target );
        if ( kmer_count > max_kmer_count )
            max_kmer_count = kmer_count;
    }

    std::string target_lca = lca.getLCA( targets );
    rep[target_lca].lca_matches++;
    read_out_lca.matches.push_back( ReadMatch{ target_lca, max_kmer_count } );
}


void classify( std::vector< Filter >&    filters,
               LCA&                      lca,
               TTax&                     tax,
               TRep&                     rep,
               SafeQueue< ReadOut >&     classified_reads_queue,
               SafeQueue< ReadOut >&     classified_lca_queue,
               SafeQueue< ReadOut >&     unclassified_reads_queue,
               Stats&                    stats,
               Config const&             config,
               SafeQueue< ReadBatches >* pointer_current,
               SafeQueue< ReadBatches >* pointer_helper,
               bool                      hierarchy_first,
               bool                      hierarchy_last,
               int16_t                   max_error_unique )
{

    // k-mer sizes should be the same among filters
    uint16_t kmer_size = filters[0].ibf.kmerSize;

    while ( true )
    {
        // Wait here until reads are available or push is over and queue is empty
        ReadBatches rb = pointer_current->pop();

        // If batch is empty exit thread
        if ( !seqan::length( rb.ids ) )
            break;

        // store unclassified reads for next iteration
        ReadBatches left_over_reads{ rb.paired };

        for ( uint32_t readID = 0; readID < seqan::length( rb.ids ); ++readID )
        {
            uint16_t read1_len = seqan::length( rb.seqs[readID] );
            // count lens just once
            if ( hierarchy_first )
                stats.sumReadLen += read1_len;

            TMatches matches;
            ReadOut  read_out( rb.ids[readID] );
            ReadOut  read_out_lca( rb.ids[readID] );
            uint32_t count_filtered_matches = 0;

            if ( read1_len >= kmer_size ) // just skip classification, add read to left over (dbs can have different
                                          // kmer sizes) or unclassified
            {
                if ( rb.paired ) // paired-end mode
                {
                    uint16_t read2_len = seqan::length( rb.seqs2[readID] );
                    if ( read2_len >= kmer_size )
                    {
                        // effective length of the pair for error calculation
                        uint16_t effective_read_len = read1_len + read2_len + 1 - kmer_size;

                        // count lens just once
                        if ( hierarchy_first )
                            stats.sumReadLen += read2_len;

                        uint16_t max_kmer_count_read = find_matches_paired(
                            matches, filters, rb.seqs[readID], rb.seqs2[readID], effective_read_len );
                        count_filtered_matches = filter_matches(
                            read_out, matches, rep, effective_read_len, max_kmer_count_read, kmer_size, config.offset );

                        if ( count_filtered_matches == 1 )
                        {
                            bool real_unique = ( max_error_unique >= 0 ) ? check_unique( read_out,
                                                                                         effective_read_len,
                                                                                         kmer_size,
                                                                                         max_error_unique,
                                                                                         config.offset,
                                                                                         tax )
                                                                         : true;
                            if ( real_unique )
                                rep[read_out.matches[0].target].unique_matches++;
                            else
                                rep[read_out.matches[0].target].lca_matches++;

                            read_out_lca = read_out; // just one match
                        }
                        else
                        {
                            lca_matches( read_out, read_out_lca, lca, rep );
                        }
                    }
                }
                else // single-end mode
                {
                    uint16_t max_kmer_count_read = find_matches( matches, filters, rb.seqs[readID] );
                    count_filtered_matches       = filter_matches(
                        read_out, matches, rep, read1_len, max_kmer_count_read, kmer_size, config.offset );

                    if ( count_filtered_matches == 1 )
                    {
                        bool real_unique =
                            ( max_error_unique >= 0 )
                                ? check_unique( read_out, read1_len, kmer_size, max_error_unique, config.offset, tax )
                                : true;
                        if ( real_unique )
                            rep[read_out.matches[0].target].unique_matches++;
                        else
                            rep[read_out.matches[0].target].lca_matches++;

                        read_out_lca = read_out; // just one match
                    }
                    else
                    {
                        lca_matches( read_out, read_out_lca, lca, rep );
                    }
                }
            }

            // If there are matches, add to printing queue
            if ( count_filtered_matches > 0 )
            {
                stats.classifiedReads += 1;
                stats.matches += count_filtered_matches;
                classified_lca_queue.push( read_out_lca );
                if ( config.output_all )
                    classified_reads_queue.push( read_out );
            }
            else if ( !hierarchy_last ) // if there is more levels, store read
            {
                seqan::appendValue( left_over_reads.ids, rb.ids[readID] );
                seqan::appendValue( left_over_reads.seqs, rb.seqs[readID] );
                if ( rb.paired )
                    seqan::appendValue( left_over_reads.seqs2, rb.seqs2[readID] );
            }
            else if ( config.output_unclassified ) // no more levels and no classification, add to
                                                   // unclassified printing queue
            {
                unclassified_reads_queue.push( read_out );
            }
        }

        // if there are more levels to classify and something was left, keep reads in memory
        if ( !hierarchy_last && seqan::length( left_over_reads.ids ) > 0 )
            pointer_helper->push( left_over_reads );
    }
}

void write_report( TRep& rep, TTax& tax, std::string output_file_rep )
{
    std::ofstream out_rep( output_file_rep );
    for ( auto const & [ target, report ] : rep )
    {
        out_rep << target << '\t' << tax[target].rank << '\t' << report.direct_matches << '\t' << report.lca_matches
                << '\t' << report.unique_matches << '\t' << tax[target].name << '\n';
    }
    out_rep.close();
}

void load_filters( std::vector< Filter >& filters, std::string hierarchy_label, Config& config )
{
    uint16_t first_kmer_size = 0;
    for ( auto const& filter_config : config.parsed_hierarchy[hierarchy_label].filters )
    {

        TIbf filter;
        TMap map;
        TTax tax;

        // ibf file
        seqan::retrieve( filter, seqan::toCString( filter_config.ibf_file ) );
        // set offset to user-defined value (1==no offset)
        filter.offset = config.offset;

        if ( first_kmer_size == 0 )
        {
            first_kmer_size = seqan::getKmerSize( filter );
        }
        else if ( first_kmer_size != seqan::getKmerSize( filter ) )
        {
            std::cerr
                << "ERROR: filters on the same hierarchy should have same k-mer size configuration. Ignoring filter: "
                << filter_config.ibf_file << std::endl;
            continue;
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
        }
        infile.close();

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
}

void print_time( GanonClassify::Config& config,
                 const StopClock&       timeGanon,
                 const StopClock&       timeLoadReads,
                 const StopClock&       timeLoadFilters,
                 const StopClock&       timeClass,
                 const StopClock&       timePrintClass,
                 const StopClock&       timePrintUnclass )
{
    using ::operator<<;

    std::cerr << "ganon-classify start time: " << timeGanon.begin() << std::endl;
    std::cerr << "Loading reads  start time: " << timeLoadReads.begin() << std::endl;
    std::cerr << "Class./ Print. start time: " << timeClass.begin() << std::endl;
    std::cerr << "Loading reads    end time: " << timeLoadReads.end() << std::endl;
    std::cerr << "Classifying      end time: " << timeClass.end() << std::endl;
    std::cerr << "Printing clas.   end time: " << timePrintClass.end() << std::endl;
    if ( config.output_unclassified )
        std::cerr << "Printing unclas. end time: " << timePrintUnclass.end() << std::endl;
    std::cerr << "ganon-classify   end time: " << timeGanon.end() << std::endl;
    std::cerr << std::endl;
    std::cerr << " - loading filters: " << timeLoadFilters.elapsed() << std::endl;
    std::cerr << " - classifying (" << config.threads_classify << "t): " << timeClass.elapsed() << std::endl;
    std::cerr << " - printing: " << timePrintClass.elapsed() << std::endl;
    if ( config.output_unclassified )
        std::cerr << " - printing unclassified" << timePrintUnclass.elapsed() << std::endl;
    std::cerr << " - total: " << timeGanon.elapsed() << std::endl;
    std::cerr << std::endl;
}

void print_stats( Stats& stats, const StopClock& timeClass )
{
    const double elapsed_classification = timeClass.elapsed();
    std::cerr << "ganon-classify processed " << stats.totalReads << " sequences (" << stats.sumReadLen / 1000000.0
              << " Mbp) in " << elapsed_classification << " seconds ("
              << ( stats.totalReads / 1000.0 ) / ( elapsed_classification / 60.0 ) << " Kseq/m, "
              << ( stats.sumReadLen / 1000000.0 ) / ( elapsed_classification / 60.0 ) << " Mbp/m)" << std::endl;
    std::cerr << " - " << stats.classifiedReads << " sequences classified ("
              << ( stats.classifiedReads / static_cast< double >( stats.totalReads ) ) * 100 << "%)" << std::endl;
    std::cerr << " - " << stats.totalReads - stats.classifiedReads << " sequences unclassified ("
              << ( ( stats.totalReads - stats.classifiedReads ) / static_cast< double >( stats.totalReads ) ) * 100
              << "%)" << std::endl;
    std::cerr << " - " << stats.matches << " matches (avg. "
              << ( stats.matches / static_cast< double >( stats.classifiedReads ) ) << " match/read)" << std::endl;
}

uint32_t parse_single_reads( seqan::SeqFileIn&                 seqFileIn,
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
}

void parse_reads( SafeQueue< detail::ReadBatches >& queue1,
                  StopClock&                        timeLoadReads,
                  Stats&                            stats,
                  Config const&                     config )
{
    for ( auto const& reads_file : config.single_reads )
    {
        seqan::SeqFileIn seqFileIn;
        if ( !seqan::open( seqFileIn, seqan::toCString( reads_file ) ) )
        {
            std::cerr << "ERROR: Unable to open the file: " << reads_file << std::endl;
            continue;
        }
        uint32_t pos = 0;
        while ( !seqan::atEnd( seqFileIn ) )
        {
            pos = seqan::position( seqFileIn );
            seqan::StringSet< seqan::CharString > ids;
            seqan::StringSet< seqan::CharString > seqs;
            try
            {
                seqan::readRecords( ids, seqs, seqFileIn, config.n_reads );
            }
            catch ( seqan::Exception const& e )
            {
                // Error occured, faulty fastq, continue to parse reads one by one from the last valid position
                std::cerr << "ERROR: Problems while reading the file in batches: " << reads_file << " [" << e.what()
                          << "]. Switching to single line parsing." << std::endl;
                stats.totalReads += parse_single_reads( seqFileIn, pos, config.n_reads, queue1 );
                break;
            }
            stats.totalReads += seqan::length( ids );
            queue1.push( detail::ReadBatches{ false, ids, seqs } );
        }
        seqan::close( seqFileIn );
    }
    // paired-reads
    if ( config.paired_reads.size() > 0 )
    {
        for ( uint16_t pair_cnt = 0; pair_cnt < config.paired_reads.size(); pair_cnt += 2 )
        {
            seqan::SeqFileIn seqFileIn1;
            seqan::SeqFileIn seqFileIn2;
            if ( !seqan::open( seqFileIn1, seqan::toCString( config.paired_reads[pair_cnt] ) ) )
            {
                std::cerr << "ERROR: Unable to open the file: " << config.paired_reads[pair_cnt] << std::endl;
                continue;
            }
            if ( !seqan::open( seqFileIn2, seqan::toCString( config.paired_reads[pair_cnt + 1] ) ) )
            {
                std::cerr << "ERROR: Unable to open the file: " << config.paired_reads[pair_cnt + 1] << std::endl;
                continue;
            }
            while ( !seqan::atEnd( seqFileIn1 ) )
            {
                seqan::StringSet< seqan::CharString > ids1;
                seqan::StringSet< seqan::CharString > seqs1;
                seqan::StringSet< seqan::CharString > ids2;
                seqan::StringSet< seqan::CharString > seqs2;
                try
                {
                    seqan::readRecords( ids1, seqs1, seqFileIn1, config.n_reads );
                    seqan::readRecords( ids2, seqs2, seqFileIn2, config.n_reads );
                }
                catch ( seqan::Exception const& e )
                {
                    std::cerr << "ERROR: " << e.what() << std::endl;
                    continue;
                }
                if ( seqan::length( ids1 ) != seqan::length( ids2 ) )
                {
                    std::cerr << "ERROR: Paired-read files do not match: " << config.paired_reads[pair_cnt] << ","
                              << config.paired_reads[pair_cnt + 1] << std::endl;
                    break;
                }
                stats.totalReads += seqan::length( ids1 );
                queue1.push( detail::ReadBatches{ true, ids1, seqs1, seqs2 } );
            }
            seqan::close( seqFileIn1 );
            seqan::close( seqFileIn2 );
        }
    }
    queue1.notify_push_over();
    timeLoadReads.stop();
}

void write_classified( SafeQueue< detail::ReadOut >& classified_reads_queue,
                       std::ofstream&                out,
                       StopClock&                    timePrintClass )
{
    while ( true )
    {
        detail::ReadOut ro = classified_reads_queue.pop();
        if ( ro.readID != "" )
        {
            for ( uint32_t i = 0; i < ro.matches.size(); ++i )
            {
                out << ro.readID << '\t' << ro.matches[i].target << '\t' << ro.matches[i].kmer_count << '\n';
            }
        }
        else
        {
            timePrintClass.stop();
            break;
        }
    }
}

void write_unclassified( SafeQueue< detail::ReadOut >& unclassified_reads_queue,
                         std::string                   out_unclassified_file,
                         StopClock&                    timePrintUnclass )
{
    std::ofstream out_unclassified( out_unclassified_file );
    while ( true )
    {
        detail::ReadOut rou = unclassified_reads_queue.pop();
        if ( rou.readID != "" )
        {
            out_unclassified << rou.readID << '\n';
        }
        else
        {
            timePrintUnclass.stop();
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

LCA pre_process_lca( TTax& tax )
{
    LCA lca;
    for ( auto const & [ target, node ] : tax )
    {
        lca.addEdge( node.parent, target );
    }
    lca.doEulerWalk();
    return lca;
}


} // namespace detail

bool run( Config config )
{
    if ( !config.validate() )
        return false;

    // Time control
    StopClock timeGanon;
    timeGanon.start();
    StopClock timeLoadReads;
    StopClock timeLoadFilters;
    StopClock timeClass;
    StopClock timePrintClass;
    StopClock timePrintClassLCA;
    StopClock timePrintUnclass;

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
    timeLoadReads.start();
    std::future< void > read_task = std::async( std::launch::async,
                                                detail::parse_reads,
                                                std::ref( queue1 ),
                                                std::ref( timeLoadReads ),
                                                std::ref( stats ),
                                                std::ref( config ) );

    // SafeQueue for printing unclassified
    SafeQueue< detail::ReadOut > unclassified_reads_queue;
    // Thread for printing unclassified reads
    std::future< void > write_unclassified_task;
    if ( config.output_unclassified && !config.output_prefix.empty() )
    {
        timePrintUnclass.start();
        write_unclassified_task = std::async( std::launch::async,
                                              detail::write_unclassified,
                                              std::ref( unclassified_reads_queue ),
                                              config.output_prefix + ".unc",
                                              std::ref( timePrintUnclass ) );
    }


    // Pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_current = &queue1; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_helper  = &queue2; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_extra;             // pointer to the queues

    uint16_t hierarchy_id   = 0;
    uint16_t hierarchy_size = config.parsed_hierarchy.size();


    // For every hierarchy level
    for ( auto const & [ hierarchy_label, hierarchy_config ] : config.parsed_hierarchy )
    {
        ++hierarchy_id;
        bool hierarchy_first = ( hierarchy_id == 1 );
        bool hierarchy_last  = ( hierarchy_id == hierarchy_size );

        timeLoadFilters.start();
        std::vector< detail::Filter > filters;
        detail::load_filters( filters, hierarchy_label, config );
        timeLoadFilters.stop();

        // Reports
        detail::TRep rep;

        // merge repeated elements
        detail::TTax tax = detail::merge_tax( filters );

        // for ( auto const& [target, nodes] : tax )
        //    std::cerr << target << "\t" << nodes.parent << "\t" << nodes.rank << "\t" << nodes.name << std::endl;

        // pre-processing of nodes
        LCA lca = detail::pre_process_lca( tax );

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
                queue1.set_max_size( -1 );
        }

        SafeQueue< detail::ReadOut > classified_reads_queue;
        SafeQueue< detail::ReadOut > classified_lca_queue;

        if ( !config.output_prefix.empty() )
        {
            if ( hierarchy_first || !config.output_hierarchy_single )
                out.open( hierarchy_config.output_file_lca );
            else // append if not first and output_hierarchy_single
                out.open( hierarchy_config.output_file_lca, std::ofstream::app );

            if ( config.output_all )
            {
                if ( hierarchy_first || !config.output_hierarchy_single )
                    out_all.open( hierarchy_config.output_file_all );
                else // append if not first and output_hierarchy_single
                    out_all.open( hierarchy_config.output_file_all, std::ofstream::app );

                // Start writing thread for all matches
                timePrintClass.start();
                write_tasks.emplace_back( std::async( std::launch::async,
                                                      detail::write_classified,
                                                      std::ref( classified_reads_queue ),
                                                      std::ref( out_all ),
                                                      std::ref( timePrintClass ) ) );
            }
        }

        // Start writing thread for lca matches
        timePrintClassLCA.start();
        write_tasks.emplace_back( std::async( std::launch::async,
                                              detail::write_classified,
                                              std::ref( classified_lca_queue ),
                                              std::ref( out ),
                                              std::ref( timePrintClassLCA ) ) );

        std::vector< std::future< void > > tasks;
        // Threads for classification
        timeClass.start();
        for ( uint16_t taskNo = 0; taskNo < config.threads_classify; ++taskNo )
        {
            tasks.emplace_back( std::async( std::launch::async,
                                            detail::classify,
                                            std::ref( filters ),
                                            std::ref( lca ),
                                            std::ref( tax ),
                                            std::ref( rep ),
                                            std::ref( classified_reads_queue ),
                                            std::ref( classified_lca_queue ),
                                            std::ref( unclassified_reads_queue ),
                                            std::ref( stats ),
                                            std::ref( config ),
                                            pointer_current,
                                            pointer_helper,
                                            hierarchy_first,
                                            hierarchy_last,
                                            hierarchy_config.max_error_unique ) );
        }

        // Wait here until classification is over
        for ( auto&& task : tasks )
        {
            task.get();
        }
        // Inform that no more reads are going to be pushed
        classified_reads_queue.notify_push_over();
        classified_lca_queue.notify_push_over();

        // Write report
        if ( !config.output_prefix.empty() )
            detail::write_report( rep, tax, hierarchy_config.output_file_rep );

        // Wait here until all files are written
        for ( auto&& task : write_tasks )
        {
            task.get();
        }

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

        timeClass.stop();
    }

    // Wait here until all unclassified reads are written
    if ( config.output_unclassified )
    {
        unclassified_reads_queue.notify_push_over();
        write_unclassified_task.get();
    }

    timeGanon.stop();

    if ( !config.quiet )
    {
        std::cerr << std::endl;
        if ( config.verbose )
        {
            detail::print_time(
                config, timeGanon, timeLoadReads, timeLoadFilters, timeClass, timePrintClass, timePrintUnclass );
        }
        detail::print_stats( stats, timeClass );
    }
    return true;
}

} // namespace GanonClassify
