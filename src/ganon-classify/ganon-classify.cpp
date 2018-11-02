#include <atomic>
#include <chrono>
#include <ctime>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <seqan/kmer.h>
#include <utils/safequeue.hpp>
#include <vector>

#include "Arguments.hpp"

inline uint16_t kmer_threshold( uint16_t readLen, uint16_t kmerSize, uint16_t max_error )
{
    uint16_t threshold = 0;
    if ( readLen > kmerSize * ( 1 + max_error ) )
        threshold = readLen - kmerSize + 1 - ( max_error * kmerSize );
    return threshold;
}

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
    seqan::CharString        readID;
    std::vector< ReadMatch > matches;
};


int main( int argc, char** argv )
{
    std::ios_base::sync_with_stdio( false ); // speed up output to STDOUT (allows buffering) ->
                                             // https://en.cppreference.com/w/cpp/io/ios_base/sync_with_stdio


    auto ganon_start = std::chrono::high_resolution_clock::now();

    // Parse arguments
    Arguments args( argc, argv );
    if ( !args.parse() )
        return 0;

    // if ( args.verbose )
    args.print();
    //////////////////////////////


    // Set output
    std::ofstream out;
    std::ofstream out_unclassified;
    if ( !args.output_file.empty() )
    { // output to a file
        out.open( args.output_file );
    }
    else
    {
        out.copyfmt( std::cout ); // STDOUT
        out.clear( std::cout.rdstate() );
        out.basic_ios< char >::rdbuf( std::cout.rdbuf() );
    }

    typedef seqan::ModifiedString<
                seqan::ModifiedString< seqan::Dna5String, seqan::ModComplementDna >,
                seqan::ModReverse
            > reversedRead;



    std::vector< std::future< void > > read_write;
    std::atomic< uint64_t >            sumReadLen      = 0;
    std::atomic< uint64_t >            classifiedReads = 0;
    uint64_t                           totalReads      = 0;

    SafeQueue< ReadBatches > queue1;
    SafeQueue< ReadBatches > queue2;

    SafeQueue< ReadOut > classified_reads_queue;
    SafeQueue< ReadOut > unclassified_reads_queue;
    bool                 finished_read = false;
    bool                 finished_clas = false;

    typedef seqan::KmerFilter< seqan::Dna5, seqan::InterleavedBloomFilter, seqan::Uncompressed > filterType;
    struct Filter
    {
        filterType                        bloom_filter;
        std::map< uint32_t, std::string > group_bin;
        uint32_t                          numberOfBins;
        uint16_t                          kmerSize;
    };

    std::chrono::duration< double >                               loading_filter_elapsed;
    std::chrono::duration< double >                               classifying_elapsed;
    std::chrono::duration< double >                               select_elapsed;
    std::chrono::duration< double >                               filter_elapsed;
    std::chrono::time_point< std::chrono::high_resolution_clock > loading_reads_start;
    std::chrono::time_point< std::chrono::high_resolution_clock > loading_reads_end;
    std::chrono::time_point< std::chrono::high_resolution_clock > general_start;
    std::chrono::time_point< std::chrono::high_resolution_clock > classifying_end;
    std::chrono::time_point< std::chrono::high_resolution_clock > printing_classified_end;
    std::chrono::time_point< std::chrono::high_resolution_clock > printing_unclassified_end;


    loading_reads_start = std::chrono::high_resolution_clock::now();
    // extra thread for reading the reads in batches
    int num_of_batches         = 1000;
    int num_of_reads_per_batch = 400;
    // num_of_batches*num_of_reads_per_batch = max. amount of reads in memory
    read_write.emplace_back(
        std::async( std::launch::async, [=, &queue1, &finished_read, &loading_reads_end, &totalReads] {
            for ( auto const& reads_file : args.reads )
            {
                seqan::SeqFileIn seqFileIn;
                if ( !seqan::open( seqFileIn, seqan::toCString( reads_file ) ) )
                {
                    std::cerr << "Unable to open " << reads_file << std::endl;
                    continue;
                }
                while ( !seqan::atEnd( seqFileIn ) )
                {
                    // std::cerr << queue1->size() << std::endl;
                    while ( queue1.size() > num_of_batches )
                    {
                        ; // spin
                    }
                    seqan::StringSet< seqan::CharString > ids;
                    seqan::StringSet< seqan::Dna5String > seqs;
                    seqan::readRecords( ids, seqs, seqFileIn, num_of_reads_per_batch );
                    totalReads += seqan::length( ids );
                    queue1.push( ReadBatches{ ids, seqs } );
                }
                seqan::close( seqFileIn );
            }
            finished_read     = true;
            loading_reads_end = std::chrono::high_resolution_clock::now();
        } ) );

    general_start = std::chrono::high_resolution_clock::now();
    // extra thread for printing classified reads
    read_write.emplace_back(
        std::async( std::launch::async, 
                    [=, &classified_reads_queue, &out, &finished_clas, &printing_classified_end] {
                        while ( true )
                        {
                            ReadOut ro = classified_reads_queue.pop();
                            for ( uint32_t i = 0; i < ro.matches.size(); ++i )
                                out << ro.readID << '\t' << ro.matches[i].group << '\t' << ro.matches[i].kmer_count << '\n';
                            if ( finished_clas && classified_reads_queue.empty() )
                            {
                                printing_classified_end = std::chrono::high_resolution_clock::now();
                                break;
                            }
                        }
                    } 
                ) 
    );

    // extra thread for printing unclassified reads
    if ( args.output_unclassified )
    {
        out_unclassified.open( args.output_unclassified_file );
        read_write.emplace_back(
            std::async( std::launch::async,
                        [=, &unclassified_reads_queue, &out_unclassified, &finished_clas, &printing_unclassified_end] {
                            while ( true )
                            {
                                ReadOut rou = unclassified_reads_queue.pop();
                                if ( rou.readID != "" ) // if not empty
                                    out_unclassified << rou.readID << '\n';
                                if ( finished_clas && unclassified_reads_queue.empty() )
                                {
                                    printing_unclassified_end = std::chrono::high_resolution_clock::now();
                                    break;
                                }
                            }
                        } ) );
    }


    SafeQueue< ReadBatches >* pointer_current; // pointer to the queues
    SafeQueue< ReadBatches >* pointer_helper;  // pointer to the queues
    SafeQueue< ReadBatches >* pointer_extra;   // pointer to the queues

    uint16_t hierarchy_id   = 0;
    uint16_t hierarchy_size = args.filters.size();
    for ( auto const& hierarchy : args.filters )
    { // filter[h] = vector<(filter, map)> -> map already sort by key
        ++hierarchy_id;
        // std::string hierarchy_name = hierarchy.first;

        std::vector< Filter > filter_hierarchy;

        auto loading_filter_start = std::chrono::high_resolution_clock::now();
        for ( auto const& file : hierarchy.second )
        {
            std::string bloom_filter_file_hierarchy = std::get< 0 >( file );
            std::string group_bin_file_hierarchy    = std::get< 1 >( file );

            // group bin files
            std::map< uint32_t, std::string > group_bin;
            std::ifstream                     infile( group_bin_file_hierarchy );
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
            filterType filter;
            seqan::retrieve( filter, seqan::toCString( bloom_filter_file_hierarchy ) );
            filter_hierarchy.push_back( Filter{
                std::move( filter ), group_bin, seqan::getNumberOfBins( filter ), seqan::getKmerSize( filter ) } );
        }
        loading_filter_elapsed += std::chrono::high_resolution_clock::now() - loading_filter_start;

        auto classifying_start = std::chrono::high_resolution_clock::now();

        // manage queues instance pointers
        if ( hierarchy_id == 1 )
        {
            pointer_current = &queue1;
            pointer_helper  = &queue2;
        }
        else
        {
            pointer_extra   = pointer_current;
            pointer_current = pointer_helper;
            pointer_helper  = pointer_extra;
        }

        // std::cerr << hierarchy_id << " - queue1 address: " << &queue1 << std::endl;
        // std::cerr << hierarchy_id << " - queue2 address: " << &queue2 << std::endl;
        // std::cerr << hierarchy_id << " - pointer_current: " << pointer_current << " - pointer_helper: " <<
        // pointer_helper << std::endl;

        std::vector< std::future< void > > tasks;

        for ( uint16_t taskNo = 0; taskNo < args.clas_threads; ++taskNo )
        {
            tasks.emplace_back( std::async( std::launch::async, [=, &filter_hierarchy,&classified_reads_queue,&unclassified_reads_queue,&finished_read,&select_elapsed,&filter_elapsed,&sumReadLen,&classifiedReads] 
            {
                while ( true )
                {
                    // std::cerr << pointer_current.size() << std::endl; //check if queue is getting empty (print 0's)
                    ReadBatches rb = pointer_current->pop();
                    if ( rb.ids != "" )
                    {                                // if not empty
                        ReadBatches left_over_reads; // store unclassified reads for next iteration
                        for ( uint32_t readID = 0; readID < seqan::length( rb.ids ); ++readID )
                        {
                            uint16_t readLen = seqan::length( rb.seqs[readID] );
                            // count lens just once
                            if ( hierarchy_id == 1 )
                                sumReadLen += readLen;

                            // k-mer sizes should be the same among filters, groups should not overlap
                            uint16_t kmerSize = filter_hierarchy[0].kmerSize;

                            std::unordered_map< std::string, int16_t > groups;
                            uint16_t                                   maxKmerCountRead = 0;
                            uint16_t threshold = kmer_threshold( readLen, kmerSize, args.max_error );

                            std::chrono::time_point< std::chrono::high_resolution_clock > filter_start;
                            // for every filter in this level
                            for ( Filter& filter : filter_hierarchy )
                            {

                                auto                    select_start = std::chrono::high_resolution_clock::now();
                                std::vector< uint32_t > selectedBins( filter.numberOfBins, 0 );
                                filter.bloom_filter.select( selectedBins, rb.seqs[readID] );
                                std::vector< uint32_t > selectedBinsRev( filter.numberOfBins, 0 );
                                filter.bloom_filter.select( selectedBinsRev, reversedRead( rb.seqs[readID] ) );
                                select_elapsed += std::chrono::high_resolution_clock::now() - select_start;

                                filter_start = std::chrono::high_resolution_clock::now();
                                for ( uint32_t binNo = 0; binNo < filter.numberOfBins; ++binNo )
                                {
                                    if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
                                    {
                                        // get best matching strand
                                        uint16_t maxKmerCountBin =
                                            std::max( selectedBins[binNo], selectedBinsRev[binNo] );
                                        // keep only the best match group/read when same groups are split in several
                                        // bins
                                        if ( groups.count( filter.group_bin[binNo] ) == 0
                                             || maxKmerCountBin > groups[filter.group_bin[binNo]] )
                                        {
                                            groups[filter.group_bin[binNo]] = maxKmerCountBin;
                                            if ( maxKmerCountBin > maxKmerCountRead )
                                                maxKmerCountRead =
                                                    maxKmerCountBin; // keep track of the max kmer count for this read
                                        }
                                    }
                                }
                                filter_elapsed += std::chrono::high_resolution_clock::now() - filter_start;
                            }

                            filter_start = std::chrono::high_resolution_clock::now();
                            // get maximum possible number of error for this read
                            // (-kmerSize+readLen-maxKmerCountRead+1)/kmerSize in a ceil formula (x + y - 1) / y
                            uint16_t max_errorRead =
                                ( ( -kmerSize + readLen - maxKmerCountRead + 1 ) + kmerSize - 1 ) / kmerSize;
                            // get min kmer count necesary to achieve the calculated number of errors
                            uint16_t minKmerCount = readLen - kmerSize + 1 - ( max_errorRead * kmerSize );

                            ReadOut classified_read_out;
                            bool    classified = false;
                            for ( auto const& v : groups )
                            { // groups[group] = kmerCount
                                if ( v.second >= minKmerCount )
                                { // apply strata filter
                                    classified_read_out.matches.push_back( ReadMatch{ v.first, v.second } );
                                    classified = true;
                                }
                            }

                            if ( classified )
                            {
                                classifiedReads += 1;
                                classified_read_out.readID = rb.ids[readID];

                                // set as negative for unique filtering
                                if ( args.unique_filtering )
                                    // if there's only one match and kmer count is lower than expected
                                    if ( classified_read_out.matches.size() == 1
                                         && classified_read_out.matches[0].kmer_count
                                                < readLen - kmerSize + 1 - ( args.max_error_unique * kmerSize ) )
                                        classified_read_out.matches[0].kmer_count =
                                            -classified_read_out.matches[0].kmer_count;

                                classified_reads_queue.push( classified_read_out );

                                // if there is more level for classification
                            }
                            else if ( hierarchy_id < hierarchy_size )
                            {
                                seqan::appendValue( left_over_reads.ids, rb.ids[readID] );
                                seqan::appendValue( left_over_reads.seqs, rb.seqs[readID] );
                            }
                            else if ( args.output_unclassified )
                            {
                                ReadOut unclassified_read_out;
                                unclassified_read_out.readID = rb.ids[readID];
                                unclassified_reads_queue.push( unclassified_read_out );
                            }
                            filter_elapsed += std::chrono::high_resolution_clock::now() - filter_start;
                        }

                        // if something was added to the classified reads (there are more levels, keep reads in memory)
                        if ( seqan::length( left_over_reads.ids ) > 0 )
                            pointer_helper->push( left_over_reads );
                    }

                    if ( finished_read && pointer_current->empty() )
                    { // if finished reading from file (first iter) and current queue is empty
                        break;
                    }
                }
            } ) );
        }
        for ( auto&& task : tasks )
        {
            task.get();
        }
        classifying_elapsed += std::chrono::high_resolution_clock::now() - classifying_start;
    }
    finished_clas   = true;
    classifying_end = std::chrono::high_resolution_clock::now();

    for ( auto&& task : read_write )
    {
        task.get();
    }

    if ( args.output_file != "" )
        out.close();

    if ( args.output_unclassified )
        out_unclassified.close();

    auto ganon_end      = std::chrono::high_resolution_clock::now();
    auto ganon_end_time = std::chrono::system_clock::to_time_t( ganon_end );

    std::cerr << std::endl;
    auto ganon_start_time = std::chrono::system_clock::to_time_t( ganon_start );
    std::cerr << "ganon-classify start time: " << std::ctime( &ganon_start_time );
    auto loading_reads_start_time = std::chrono::system_clock::to_time_t( loading_reads_start );
    std::cerr << "Loading reads  start time: " << std::ctime( &loading_reads_start_time );
    auto general_start_time = std::chrono::system_clock::to_time_t( general_start );
    std::cerr << "Class./ Print. start time: " << std::ctime( &general_start_time );
    auto loading_reads_end_time = std::chrono::system_clock::to_time_t( loading_reads_end );
    std::cerr << "Loading reads    end time: " << std::ctime( &loading_reads_end_time );
    auto classifying_end_time = std::chrono::system_clock::to_time_t( classifying_end );
    std::cerr << "Classifying      end time: " << std::ctime( &classifying_end_time );
    auto printing_classified_end_time = std::chrono::system_clock::to_time_t( printing_classified_end );
    std::cerr << "Printing clas.   end time: " << std::ctime( &printing_classified_end_time );
    if ( args.output_unclassified )
    {
        auto printing_unclassified_end_time = std::chrono::system_clock::to_time_t( printing_unclassified_end );
        std::cerr << "Printing unclas. end time: " << std::ctime( &printing_unclassified_end_time );
    }
    std::cerr << "ganon-classify   end time: " << std::ctime( &ganon_end_time ) << std::endl;


    std::cerr << " 1) loading filters: " << loading_filter_elapsed.count() << std::endl;
    double total_classifying_elapsed = classifying_elapsed.count();
    std::cerr << " 2) classifying (" << args.clas_threads << "t): " << total_classifying_elapsed << std::endl;
    std::cerr << "    - select: "
              << total_classifying_elapsed
                     * ( select_elapsed.count() / ( select_elapsed.count() + filter_elapsed.count() ) )
              << " (" << select_elapsed.count() << ")" << std::endl;
    std::cerr << "    - filter: "
              << total_classifying_elapsed
                     * ( filter_elapsed.count() / ( select_elapsed.count() + filter_elapsed.count() ) )
              << " (" << filter_elapsed.count() << ")" << std::endl;
    std::cerr << " Total elapsed: " << std::chrono::duration< double >( ganon_end - ganon_start ).count() << std::endl;
    std::cerr << std::endl;

    std::cerr << "ganon-classify processed " << totalReads << " sequences (" << sumReadLen / 1000000.0 << " Mbp) in "
              << total_classifying_elapsed << " seconds ("
              << ( totalReads / 1000.0 ) / ( total_classifying_elapsed / 60.0 ) << " Kseq/m, "
              << ( sumReadLen / 1000000.0 ) / ( total_classifying_elapsed / 60.0 ) << " Mbp/m)" << std::endl;
    std::cerr << " - " << classifiedReads << " sequences classified ("
              << ( classifiedReads / (double) totalReads ) * 100 << "%)" << std::endl;
    std::cerr << " - " << totalReads - classifiedReads << " sequences unclassified ("
              << ( ( totalReads - classifiedReads ) / (double) totalReads ) * 100 << "%)" << std::endl;

    return 0;
}
