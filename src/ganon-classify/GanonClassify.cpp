#include "GanonClassify.hpp"

#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

#include <seqan/binning_directory.h>

#include <atomic>
#include <chrono>
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
    TBloomFilter;

typedef seqan::ModifiedString< seqan::ModifiedString< seqan::Dna5String, seqan::ModComplementDna >, seqan::ModReverse >
    TSeqRevComp;

typedef std::unordered_map< std::string, int16_t > TMatches;


struct ReadBatches
{
    bool                                  paired = false;
    seqan::StringSet< seqan::CharString > ids;
    seqan::StringSet< seqan::Dna5String > seqs;
    seqan::StringSet< seqan::Dna5String > seqs2;
};

struct ReadMatch
{
    std::string group;
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
    TBloomFilter                      bloom_filter;
    std::map< uint32_t, std::string > group_bin;
    uint32_t                          numberOfBins;
    uint16_t                          kmerSize;
    FilterConfig                      filter_config;
};

inline uint16_t get_error( uint16_t readLen, uint16_t kmerSize, uint16_t kmer_count, uint16_t offset )
{
    // Return the optimal number of errors for a certain sequence based on the kmer_count
    // (offset-1) -> to correct for the floor left overs
    return std::ceil( ( -kmerSize + readLen - ( kmer_count * offset + ( offset - 1 ) ) + 1 ) / (float) kmerSize );
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
    // Return threshold (number of kmers) based on an percentage of kmers
    // ceil -> round-up min # k-mers, floor -> round-down for offset
    return std::floor( std::ceil( ( readLen - kmerSize + 1 ) * min_kmers ) / offset );
}

inline void flag_max_error_unique( ReadOut& read_out, uint16_t threshold_error_unique )
{
    // set as negative for unique filtering
    // if kmer count is lower than expected
    if ( read_out.matches[0].kmer_count < threshold_error_unique )
        read_out.matches[0].kmer_count = -read_out.matches[0].kmer_count;
}

void select_matches( TMatches&                matches,
                     std::vector< uint16_t >& selectedBins,
                     std::vector< uint16_t >& selectedBinsRev,
                     Filter&                  filter,
                     uint16_t                 threshold,
                     uint16_t&                maxKmerCountRead )
{
    // for each bin
    for ( uint32_t binNo = 0; binNo < filter.numberOfBins; ++binNo )
    {
        // if kmer count is higher than threshold
        if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
        {
            // get best matching strand
            uint16_t maxKmerCountBin = std::max( selectedBins[binNo], selectedBinsRev[binNo] );
            // keep only the best match group/read when same groups are split in several
            // bins
            if ( matches.count( filter.group_bin[binNo] ) == 0 || maxKmerCountBin > matches[filter.group_bin[binNo]] )
            {
                // store match to group
                matches[filter.group_bin[binNo]] = maxKmerCountBin;
                if ( maxKmerCountBin > maxKmerCountRead )
                    maxKmerCountRead = maxKmerCountBin;
            }
        }
    }
}

uint16_t find_matches( TMatches& matches, std::vector< Filter >& filter_hierarchy, seqan::Dna5String& read_seq )
{
    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;
    for ( Filter& filter : filter_hierarchy )
    {
        // IBF count
        std::vector< uint16_t > selectedBins    = seqan::count( filter.bloom_filter, read_seq );
        std::vector< uint16_t > selectedBinsRev = seqan::count( filter.bloom_filter, TSeqRevComp( read_seq ) );

        uint16_t threshold = ( filter.filter_config.min_kmers > 0 )
                                 ? get_threshold_kmers( seqan::length( read_seq ),
                                                        filter.bloom_filter.kmerSize,
                                                        filter.filter_config.min_kmers,
                                                        filter.bloom_filter.offset )
                                 : get_threshold_errors( seqan::length( read_seq ),
                                                         filter.bloom_filter.kmerSize,
                                                         filter.filter_config.max_error,
                                                         filter.bloom_filter.offset );

        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filter, threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint16_t find_matches_paired( TMatches&              matches,
                              std::vector< Filter >& filter_hierarchy,
                              seqan::Dna5String&     read_seq,
                              seqan::Dna5String&     read_seq2,
                              uint16_t               effective_read_len )
{
    // send it as a reference to be valid for every filter
    uint16_t max_kmer_count_read = 0;
    for ( Filter& filter : filter_hierarchy )
    {
        // IBF count
        // FR
        std::vector< uint16_t > selectedBins = seqan::count( filter.bloom_filter, read_seq );
        filter.bloom_filter.count< THashCount >( selectedBins, TSeqRevComp( read_seq2 ) );
        // RF
        std::vector< uint16_t > selectedBinsRev = seqan::count( filter.bloom_filter, TSeqRevComp( read_seq ) );
        filter.bloom_filter.count< THashCount >( selectedBinsRev, read_seq2 );

        uint16_t threshold = ( filter.filter_config.min_kmers > 0 )
                                 ? get_threshold_kmers( effective_read_len,
                                                        filter.bloom_filter.kmerSize,
                                                        filter.filter_config.min_kmers,
                                                        filter.bloom_filter.offset )
                                 : get_threshold_errors( effective_read_len,
                                                         filter.bloom_filter.kmerSize,
                                                         filter.filter_config.max_error,
                                                         filter.bloom_filter.offset );

        // select matches above chosen threshold
        select_matches( matches, selectedBins, selectedBinsRev, filter, threshold, max_kmer_count_read );
    }
    return max_kmer_count_read;
}

uint32_t filter_matches( ReadOut&  read_out,
                         TMatches& matches,
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
    { // matches[group] = kmerCount
        if ( v.second >= threshold_strata )
        { // apply strata filter
            read_out.matches.push_back( ReadMatch{ v.first, v.second } );
        }
    }

    return read_out.matches.size();
}


void classify( std::vector< Filter >&    filter_hierarchy,
               SafeQueue< ReadOut >&     classified_reads_queue,
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
    uint16_t kmer_size = filter_hierarchy[0].kmerSize;

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
                            matches, filter_hierarchy, rb.seqs[readID], rb.seqs2[readID], effective_read_len );
                        count_filtered_matches = filter_matches(
                            read_out, matches, effective_read_len, max_kmer_count_read, kmer_size, config.offset );

                        if ( max_error_unique >= 0 && count_filtered_matches == 1 )
                        {
                            uint16_t threshold_error_unique =
                                get_threshold_errors( effective_read_len, kmer_size, max_error_unique, config.offset );
                            flag_max_error_unique( read_out, threshold_error_unique );
                        }
                    }
                }
                else // single-end mode
                {
                    uint16_t max_kmer_count_read = find_matches( matches, filter_hierarchy, rb.seqs[readID] );
                    count_filtered_matches =
                        filter_matches( read_out, matches, read1_len, max_kmer_count_read, kmer_size, config.offset );
                    if ( max_error_unique >= 0 && count_filtered_matches == 1 )
                    {
                        uint16_t threshold_error_unique =
                            get_threshold_errors( read1_len, kmer_size, max_error_unique, config.offset );
                        flag_max_error_unique( read_out, threshold_error_unique );
                    }
                }
            }

            // If there are matches, add to printing queue
            if ( count_filtered_matches > 0 )
            {
                stats.classifiedReads += 1;
                stats.matches += count_filtered_matches;
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

void load_filters( std::vector< Filter >& filter_hierarchy, std::string hierarchy_name, Config& config )
{
    for ( auto const& filter_config : config.h_filters[hierarchy_name].filters )
    {
        // group bin files
        std::map< uint32_t, std::string > group_bin;
        std::ifstream                     infile( filter_config.group_bin_file );
        std::string                       line;
        while ( std::getline( infile, line, '\n' ) )
        {
            std::istringstream         stream_line( line );
            std::vector< std::string > fields;
            std::string                field;
            while ( std::getline( stream_line, field, '\t' ) )
                fields.push_back( field );
            // group <tab> binid
            group_bin[std::stoul( fields[1] )] = fields[0];
        }

        // bloom filter
        TBloomFilter filter;
        seqan::retrieve( filter, seqan::toCString( filter_config.bloom_filter_file ) );
        // set offset to user-defined value (1==no offset)
        filter.offset = config.offset;

        filter_hierarchy.push_back( Filter{ std::move( filter ),
                                            group_bin,
                                            seqan::getNumberOfBins( filter ),
                                            seqan::getKmerSize( filter ),
                                            filter_config } );
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
    std::cerr << " - classifying (" << config.clas_threads << "t): " << timeClass.elapsed() << std::endl;
    std::cerr << " - printing (1t): " << timePrintClass.elapsed() << std::endl;
    if ( config.output_unclassified )
        std::cerr << " - printing unclassified (1t) " << timePrintUnclass.elapsed() << std::endl;
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
              << ( stats.classifiedReads / (double) stats.totalReads ) * 100 << "%)" << std::endl;
    std::cerr << " - " << stats.totalReads - stats.classifiedReads << " sequences unclassified ("
              << ( ( stats.totalReads - stats.classifiedReads ) / (double) stats.totalReads ) * 100 << "%)"
              << std::endl;
    std::cerr << " - " << stats.matches << " matches (avg. " << ( stats.matches / (double) stats.classifiedReads )
              << " match/read)" << std::endl;
}

uint32_t parse_reads_single( seqan::SeqFileIn&                 seqFileIn,
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
    for ( auto const& reads_file : config.reads_single )
    {
        seqan::SeqFileIn seqFileIn;
        if ( !seqan::open( seqFileIn, seqan::toCString( reads_file ) ) )
        {
            std::cerr << "Unable to open " << reads_file << std::endl;
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
                std::cerr << "Error while reading file in batches: " << reads_file << " [" << e.what()
                          << "]. Switching to single line parsing." << std::endl;
                stats.totalReads += parse_reads_single( seqFileIn, pos, config.n_reads, queue1 );
                break;
            }
            stats.totalReads += seqan::length( ids );
            queue1.push( detail::ReadBatches{ false, ids, seqs } );
        }
        seqan::close( seqFileIn );
    }
    // paired-reads
    if ( config.reads_paired.size() > 0 )
    {
        for ( uint16_t pair_cnt = 0; pair_cnt < config.reads_paired.size(); pair_cnt += 2 )
        {
            seqan::SeqFileIn seqFileIn1;
            seqan::SeqFileIn seqFileIn2;
            if ( !seqan::open( seqFileIn1, seqan::toCString( config.reads_paired[pair_cnt] ) ) )
            {
                std::cerr << "Unable to open " << config.reads_paired[pair_cnt] << std::endl;
                continue;
            }
            if ( !seqan::open( seqFileIn2, seqan::toCString( config.reads_paired[pair_cnt + 1] ) ) )
            {
                std::cerr << "Unable to open " << config.reads_paired[pair_cnt + 1] << std::endl;
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
                    std::cerr << e.what() << std::endl;
                    continue;
                }
                if ( seqan::length( ids1 ) != seqan::length( ids2 ) )
                {
                    std::cerr << "Paired-read files do not match" << std::endl;
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


} // namespace detail

bool run( Config config )
{
    // speed up output to STDOUT (allows buffering)
    std::ios_base::sync_with_stdio( false );

    if ( !config.validate() )
        return false;

    // Time control
    StopClock timeGanon;
    timeGanon.start();
    StopClock timeLoadReads;
    StopClock timeLoadFilters;
    StopClock timeClass;
    StopClock timePrintClass;
    StopClock timePrintUnclass;

    if ( config.verbose )
    {
        std::cerr << config;
    }

    //////////////////////////////

    // Set output stream (file or stdout)
    std::ofstream out;
    std::ofstream out_unclassified;
    // If there's no output file, redirect to STDOUT
    if ( config.output_file.empty() )
    {
        out.copyfmt( std::cout ); // STDOUT
        out.clear( std::cout.rdstate() );
        out.basic_ios< char >::rdbuf( std::cout.rdbuf() );
    }
    else if ( !config.split_output_file_hierarchy )
    {
        out.open( config.output_file );
    }

    // Queues for internal read handling
    // queue1 get reads from file
    // queue2 will get unclassified reads if hierachy == 2
    // if hierachy == 3 queue1 is used for unclassified and so on
    SafeQueue< detail::ReadBatches > queue1(
        config.n_batches ); // config.n_batches*config.n_reads = max. amount of reads in memory
    SafeQueue< detail::ReadBatches > queue2;

    // Queues for classified, unclassified reads (print)
    SafeQueue< detail::ReadOut > classified_reads_queue;
    SafeQueue< detail::ReadOut > unclassified_reads_queue;

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


    // Thread for printing classified reads
    timePrintClass.start();
    std::vector< std::future< void > > write_tasks;
    write_tasks.emplace_back( std::async( std::launch::async, [=, &classified_reads_queue, &out, &timePrintClass] {
        while ( true )
        {
            detail::ReadOut ro = classified_reads_queue.pop();
            if ( ro.readID != "" )
            {
                for ( uint32_t i = 0; i < ro.matches.size(); ++i )
                {
                    out << ro.readID << '\t' << ro.matches[i].group << '\t' << ro.matches[i].kmer_count << '\n';
                }
            }
            else
            {
                timePrintClass.stop();
                break;
            }
        }
    } ) );

    // Thread for printing unclassified reads
    timePrintUnclass.start();
    if ( config.output_unclassified )
    {
        out_unclassified.open( config.output_unclassified_file );
        write_tasks.emplace_back(
            std::async( std::launch::async, [=, &unclassified_reads_queue, &out_unclassified, &timePrintUnclass] {
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
                        break;
                    }
                }
            } ) );
    }


    // Pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_current = &queue1; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_helper  = &queue2; // pointer to the queues
    SafeQueue< detail::ReadBatches >* pointer_extra;             // pointer to the queues

    uint16_t hierarchy_id   = 0;
    uint16_t hierarchy_size = config.h_filters.size();
    // For every hierarchy level
    for ( auto const& hierarchy : config.h_filters )
    {
        ++hierarchy_id;
        bool hierarchy_first = ( hierarchy_id == 1 );
        bool hierarchy_last  = ( hierarchy_id == hierarchy_size );

        std::string hierarchy_name = hierarchy.first;

        // open hierarchy output file (before filter, ganon wrapper is watching)
        if ( config.split_output_file_hierarchy && !hierarchy.second.output_file.empty() )
            out.open( hierarchy.second.output_file );

        timeLoadFilters.start();
        std::vector< detail::Filter > filter_hierarchy;
        detail::load_filters( filter_hierarchy, hierarchy_name, config );
        timeLoadFilters.stop();

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


        std::vector< std::future< void > > tasks;
        // Threads for classification
        timeClass.start();
        for ( uint16_t taskNo = 0; taskNo < config.clas_threads; ++taskNo )
        {
            tasks.emplace_back( std::async( std::launch::async,
                                            detail::classify,
                                            std::ref( filter_hierarchy ),
                                            std::ref( classified_reads_queue ),
                                            std::ref( unclassified_reads_queue ),
                                            std::ref( stats ),
                                            std::ref( config ),
                                            pointer_current,
                                            pointer_helper,
                                            hierarchy_first,
                                            hierarchy_last,
                                            hierarchy.second.max_error_unique ) );
        }


        for ( auto&& task : tasks )
        {
            task.get();
        }

        if ( hierarchy_first )
        {
            read_task.get();                    // get reading tasks at the end of the first hierarchy
            pointer_helper->notify_push_over(); // notify push is over, only on first time (will be always set over for
                                                // next iterations)
        }

        if ( config.split_output_file_hierarchy && !hierarchy.second.output_file.empty() )
        {
            while ( !classified_reads_queue.empty() )
                ; // spin
            // wait to have all writen before closing - may still fail for the last entry (wait a second) - fix next
            // version
            std::chrono::microseconds ms( 500000 ); // half second
            std::this_thread::sleep_for( ms );
            out << '\n'; // write line break at the end, signal to ganon wrapper that file is over in case of multiple
                         // files
            out.close();
        }

        timeClass.stop();
    }

    // notify that classification stoped adding new items, can exit when finished
    classified_reads_queue.notify_push_over();
    unclassified_reads_queue.notify_push_over();
    for ( auto&& task : write_tasks )
    {
        task.get();
    }

    if ( config.output_unclassified )
        out_unclassified.close();
    if ( !config.split_output_file_hierarchy && !config.output_file.empty() )
    {
        out << '\n'; // write line break at the end, signal to ganon wrapper that file is over in case of multiple files
        out.close();
    }

    timeGanon.stop();

    std::cerr << std::endl;
    if ( config.verbose )
    {
        detail::print_time(
            config, timeGanon, timeLoadReads, timeLoadFilters, timeClass, timePrintClass, timePrintUnclass );
    }
    detail::print_stats( stats, timeClass );

    return true;
}

} // namespace GanonClassify
