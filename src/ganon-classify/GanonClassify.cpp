#include "GanonClassify.hpp"

#include <utils/SafeQueue.hpp>
#include <utils/Time.hpp>

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

#ifdef GANON_OFFSET
// Filter is created with seqan::Offset<1>
// when using seqan::count it acts like seqan::Normal if offset=1
typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                 seqan::BDConfig< seqan::Dna5, seqan::Offset< 1 >, seqan::Uncompressed > >
    Tfilter;
#else
typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                 seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
    Tfilter;
#endif

typedef seqan::ModifiedString< seqan::ModifiedString< seqan::Dna5String, seqan::ModComplementDna >, seqan::ModReverse >
    reversedRead;

typedef std::unordered_map< std::string, int16_t > Tmatches;


struct ReadBatches
{
    seqan::StringSet< seqan::CharString > ids;
    seqan::StringSet< seqan::Dna5String > seqs;
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
    Tfilter                           bloom_filter;
    std::map< uint32_t, std::string > group_bin;
    uint32_t                          numberOfBins;
    uint16_t                          kmerSize;
    FilterConfig                      filter_config;
};

inline uint16_t get_threshold( uint16_t readLen, uint16_t kmerSize, uint16_t max_error, uint16_t offset )
{
    return readLen + 1u > kmerSize * ( 1u + max_error )
               ? std::floor( ( readLen - kmerSize + 1u - ( max_error * kmerSize ) ) / offset )
               : 0u;
}

inline uint16_t classify_read( Tmatches& matches, std::vector< Filter >& filter_hierarchy, seqan::Dna5String& read_seq )
{
    // std::chrono::time_point< std::chrono::high_resolution_clock > filter_start;
    // for every filter in this level
    uint16_t maxKmerCountRead = 0;
    for ( Filter& filter : filter_hierarchy )
    {
        // should match threshold from lib:
        // uint32_t                threshold       = filter.filter_config.max_error;
        // std::vector< uint16_t > selectedBins    = seqan::count( filter.bloom_filter, read_seq, threshold );
        uint16_t threshold = get_threshold( seqan::length( read_seq ),
                                            filter.bloom_filter.kmerSize,
                                            filter.filter_config.max_error,
                                            filter.bloom_filter.offset );

        // if threshold == 0 it is not possible to confidently classify read given kmerSize and max_errors values
        // TODO add warning to user
        if ( threshold == 0 )
            continue;

        std::vector< uint16_t > selectedBins    = seqan::count( filter.bloom_filter, read_seq );
        std::vector< uint16_t > selectedBinsRev = seqan::count( filter.bloom_filter, reversedRead( read_seq ) );

        // select_elapsed += std::chrono::high_resolution_clock::now() - select_start;

        // filter_start = std::chrono::high_resolution_clock::now();
        for ( uint32_t binNo = 0; binNo < filter.numberOfBins; ++binNo )
        {
            if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
            {
                // get best matching strand
                uint16_t maxKmerCountBin = std::max( selectedBins[binNo], selectedBinsRev[binNo] );
                // keep only the best match group/read when same groups are split in several
                // bins
                if ( matches.count( filter.group_bin[binNo] ) == 0
                     || maxKmerCountBin > matches[filter.group_bin[binNo]] )
                {
                    matches[filter.group_bin[binNo]] = maxKmerCountBin;
                    if ( maxKmerCountBin > maxKmerCountRead )
                        maxKmerCountRead = maxKmerCountBin; // keep track of the max kmer
                                                            // count for this read
                }
            }
        }
        // filter_elapsed += std::chrono::high_resolution_clock::now() - filter_start;
    }

    return maxKmerCountRead;
}

inline uint32_t filter_matches( ReadOut&                     read_out,
                                Tmatches&                    matches,
                                uint16_t                     kmerSize,
                                uint16_t                     readLen,
                                uint16_t                     maxKmerCountRead,
                                GanonClassify::Config const& config,
                                int16_t                      max_error_unique )
{


    // get maximum possible number of error for this read
    // (config.offset-1) -> to correct for the floor left overs
    uint16_t max_errorRead = std::ceil(
        ( -kmerSize + readLen - ( maxKmerCountRead * config.offset + ( config.offset - 1 ) ) + 1 ) / (float) kmerSize );

    // get min kmer count necesary to achieve the calculated number of errors
    uint16_t threshold_strata = get_threshold( readLen, kmerSize, max_errorRead, config.offset );

    for ( auto const& v : matches )
    { // matches[group] = kmerCount
        if ( v.second >= threshold_strata )
        { // apply strata filter
            read_out.matches.push_back( ReadMatch{ v.first, v.second } );
        }
    }

    // set as negative for unique filtering
    if ( max_error_unique >= 0 )
        // if there's only one match and kmer count is lower than expected
        if ( read_out.matches.size() == 1
             && read_out.matches[0].kmer_count < get_threshold( readLen, kmerSize, max_error_unique, config.offset ) )
            read_out.matches[0].kmer_count = -read_out.matches[0].kmer_count;

    return read_out.matches.size();
}

void classify( std::vector< Filter >&    filter_hierarchy,
               SafeQueue< ReadOut >&     classified_reads_queue,
               SafeQueue< ReadOut >&     unclassified_reads_queue,
               Stats&                    stats,
               Config const&             config,
               SafeQueue< ReadBatches >* pointer_current,
               SafeQueue< ReadBatches >* pointer_helper,
               uint16_t                  hierarchy_id,
               uint16_t                  hierarchy_size,
               int16_t                   max_error_unique )
{

    while ( true )
    {
        ReadBatches rb = pointer_current->pop();
        // std::cerr << "Queue size "  << pointer_current->size() << std::endl; //check if queue is getting empty (print
        // 0's)
        if ( seqan::length( rb.ids ) > 0 )
        {
            ReadBatches left_over_reads; // store unclassified reads for next iteration
            for ( uint32_t readID = 0; readID < seqan::length( rb.ids ); ++readID )
            {
                uint16_t readLen = seqan::length( rb.seqs[readID] );
                // count lens just once
                if ( hierarchy_id == 1 )
                    stats.sumReadLen += readLen;

                // k-mer sizes should be the same among filters, groups should not overlap
                uint16_t kmerSize = filter_hierarchy[0].kmerSize;

                // Classify reads, returing matches and max kmer count achieved for the read
                Tmatches matches;
                uint16_t maxKmerCountRead = classify_read( matches, filter_hierarchy, rb.seqs[readID] );

                // Filter reads by number of errors, special filtering for unique matches
                ReadOut  read_out( rb.ids[readID] );
                uint32_t count_filtered_matches =
                    filter_matches( read_out, matches, kmerSize, readLen, maxKmerCountRead, config, max_error_unique );

                // If there are matches, add to printing queue
                if ( count_filtered_matches > 0 )
                {
                    stats.classifiedReads += 1;
                    stats.matches += count_filtered_matches;
                    classified_reads_queue.push( read_out );
                }
                else if ( hierarchy_id < hierarchy_size ) // if there is more levels, store read
                {
                    seqan::appendValue( left_over_reads.ids, rb.ids[readID] );
                    seqan::appendValue( left_over_reads.seqs, rb.seqs[readID] );
                }
                else if ( config.output_unclassified ) // no more levels and no classification, add to
                                                       // unclassified printing queue
                {
                    unclassified_reads_queue.push( read_out );
                }
            }

            // if there are more levels to classify and something was left, keep reads in memory
            if ( hierarchy_id < hierarchy_size && seqan::length( left_over_reads.ids ) > 0 )
                pointer_helper->push( left_over_reads );
        }
        else
        {
            break;
        }
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
        Tfilter filter;
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
                 Time&                  timeGanon,
                 Time&                  timeLoadReads,
                 Time&                  timeLoadFilters,
                 Time&                  timeClass,
                 Time&                  timePrintClass,
                 Time&                  timePrintUnclass )
{
    std::cerr << "ganon-classify start time: " << timeGanon.get_start_ctime();
    std::cerr << "Loading reads  start time: " << timeLoadReads.get_start_ctime();
    std::cerr << "Class./ Print. start time: " << timeClass.get_start_ctime();
    std::cerr << "Loading reads    end time: " << timeLoadReads.get_end_ctime();
    std::cerr << "Classifying      end time: " << timeClass.get_end_ctime();
    std::cerr << "Printing clas.   end time: " << timePrintClass.get_end_ctime();
    if ( config.output_unclassified )
        std::cerr << "Printing unclas. end time: " << timePrintUnclass.get_end_ctime();
    std::cerr << "ganon-classify   end time: " << timeGanon.get_end_ctime();
    std::cerr << std::endl;
    std::cerr << " - loading filters: " << timeLoadFilters.get_elapsed() << std::endl;
    std::cerr << " - classifying (" << config.clas_threads << "t): " << timeClass.get_elapsed() << std::endl;
    std::cerr << " - printing (1t): " << timePrintClass.get_elapsed() << std::endl;
    if ( config.output_unclassified )
        std::cerr << " - printing unclassified (1t) " << timePrintUnclass.get_elapsed() << std::endl;
    std::cerr << " - total: " << timeGanon.get_elapsed() << std::endl;
    std::cerr << std::endl;
}

void print_stats( Stats& stats, Time& timeClass )
{
    const double elapsed_classification = timeClass.get_elapsed();
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

} // namespace detail


bool run( Config config )
{
    // speed up output to STDOUT (allows buffering)
    std::ios_base::sync_with_stdio( false );

    if ( !config.validate() )
        return false;

    // Time control
    Time timeGanon;
    timeGanon.start();
    Time timeLoadReads;
    Time timeLoadFilters;
    Time timeClass;
    Time timePrintClass;
    Time timePrintUnclass;

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
    std::future< void > read_task( std::async( std::launch::async, [=, &queue1, &timeLoadReads, &stats] {
        for ( auto const& reads_file : config.reads )
        {
            seqan::SeqFileIn seqFileIn;
            if ( !seqan::open( seqFileIn, seqan::toCString( reads_file ) ) )
            {
                std::cerr << "Unable to open " << reads_file << std::endl;
                continue;
            }
            while ( !seqan::atEnd( seqFileIn ) )
            {
                seqan::StringSet< seqan::CharString > ids;
                seqan::StringSet< seqan::Dna5String > seqs;
                seqan::readRecords( ids, seqs, seqFileIn, config.n_reads );
                stats.totalReads += seqan::length( ids );
                queue1.push( detail::ReadBatches{ ids, seqs } );
            }
            seqan::close( seqFileIn );
        }
        queue1.notify_push_over();
        timeLoadReads.end();
    } ) );

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
                timePrintClass.end();
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
                        timePrintUnclass.end();
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
        std::string hierarchy_name = hierarchy.first;

        timeLoadFilters.start();
        std::vector< detail::Filter > filter_hierarchy;
        detail::load_filters( filter_hierarchy, hierarchy_name, config );
        timeLoadFilters.end();

        if ( config.split_output_file_hierarchy && !hierarchy.second.output_file.empty() )
            out.open( hierarchy.second.output_file );

        // hierarchy_id = 1
        //  pointer_current=queue1, data comes from file in a limited size queue
        //  pointer_helper=queue2, empty
        // hierarchy_id > 1
        //  pointer_current=queue2, with all data already in from last iteration
        //  pointer_helper=queue1, empty
        // Exchange queues instance pointers for each hierachy (if not first)
        if ( hierarchy_id > 1 )
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
                                            hierarchy_id,
                                            hierarchy_size,
                                            hierarchy.second.max_error_unique ) );
        }

        for ( auto&& task : tasks )
        {
            task.get();
        }

        if ( hierarchy_id == 1 )
        {
            read_task.get();                    // get reading tasks at the end of the first hierarchy
            pointer_helper->notify_push_over(); // notify push is over, only on first time (will be always set over for
                                                // next iterations)
        }

        if ( config.split_output_file_hierarchy && !hierarchy.second.output_file.empty() )
            out.close();

        timeClass.end();
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
        out.close();

    timeGanon.end();

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
