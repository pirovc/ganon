#include "GanonClassify.hpp"

#include <utils/SafeQueue.hpp>
#include <utils/Time.hpp>

#include <seqan/kmer.h>

#include <atomic>
#include <cinttypes>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>


namespace detail
{

typedef seqan::KmerFilter< seqan::Dna5, seqan::InterleavedBloomFilter, seqan::Uncompressed > Tfilter;

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
};


inline uint16_t kmer_threshold( uint16_t readLen, uint16_t kmerSize, uint16_t max_error )
{
    uint16_t threshold = 0;
    if ( readLen > kmerSize * ( 1 + max_error ) )
        threshold = readLen - kmerSize + 1 - ( max_error * kmerSize );
    return threshold;
}

inline uint16_t classify_read( Tmatches&              matches,
                               std::vector< Filter >& filter_hierarchy,
                               seqan::Dna5String&     read_seq,
                               uint16_t               threshold )
{
    // std::chrono::time_point< std::chrono::high_resolution_clock > filter_start;
    // for every filter in this level
    uint16_t maxKmerCountRead = 0;
    for ( Filter& filter : filter_hierarchy )
    {
        // auto                    select_start = std::chrono::high_resolution_clock::now();
        std::vector< uint32_t > selectedBins( filter.numberOfBins, 0 );
        std::vector< uint32_t > selectedBinsRev( filter.numberOfBins, 0 );
        filter.bloom_filter.select( selectedBins, read_seq );
        filter.bloom_filter.select( selectedBinsRev, reversedRead( read_seq ) );
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

inline uint32_t filter_matches( ReadOut&               read_out,
                                Tmatches&              matches,
                                uint16_t               kmerSize,
                                uint16_t               readLen,
                                uint16_t               maxKmerCountRead,
                                GanonClassify::Config& config )
{
    // get maximum possible number of error for this read
    // (-kmerSize+readLen-maxKmerCountRead+1)/kmerSize in a ceil formula (x + y - 1) / y
    uint16_t max_errorRead = ( ( -kmerSize + readLen - maxKmerCountRead + 1 ) + kmerSize - 1 ) / kmerSize;
    // get min kmer count necesary to achieve the calculated number of errors
    uint16_t minKmerCount = readLen - kmerSize + 1 - ( max_errorRead * kmerSize );

    for ( auto const& v : matches )
    { // matches[group] = kmerCount
        if ( v.second >= minKmerCount )
        { // apply strata filter
            read_out.matches.push_back( ReadMatch{ v.first, v.second } );
        }
    }

    // set as negative for unique filtering
    if ( config.unique_filtering )
        // if there's only one match and kmer count is lower than expected
        if ( read_out.matches.size() == 1
             && read_out.matches[0].kmer_count < readLen - kmerSize + 1 - ( config.max_error_unique * kmerSize ) )
            read_out.matches[0].kmer_count = -read_out.matches[0].kmer_count;

    return read_out.matches.size();
}

void load_filters( std::vector< Filter >&                                filter_hierarchy,
                   std::vector< std::tuple< std::string, std::string > > hierarchy_files )
{
    for ( auto const& file : hierarchy_files )
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
        Tfilter filter;
        seqan::retrieve( filter, seqan::toCString( bloom_filter_file_hierarchy ) );
        filter_hierarchy.push_back(
            Filter{ std::move( filter ), group_bin, seqan::getNumberOfBins( filter ), seqan::getKmerSize( filter ) } );
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

std::vector< std::string > split( const std::string& s, char delimiter )
{
    std::vector< std::string > tokens;
    std::string                token;
    std::istringstream         tokenStream( s );
    while ( std::getline( tokenStream, token, delimiter ) )
        tokens.push_back( token );
    return tokens;
}

bool parse_hierarchy( GanonClassify::Config& config )
{

    if ( config.filter_hierarchy.empty() )
    {
        if ( config.bloom_filter_files.size() != config.group_bin_files.size() )
        {
            std::cerr << "Filter and group-bin files do not match" << std::endl;
            return false;
        }
        else
        {
            for ( uint16_t h = 0; h < config.bloom_filter_files.size(); ++h )
            {
                config.filters["1"].push_back(
                    std::make_tuple( config.bloom_filter_files[h], config.group_bin_files[h] ) );
            }
        }
    }
    else
    {
        std::vector< std::string > hierarchy = split( config.filter_hierarchy, ',' );
        if ( hierarchy.size() != config.bloom_filter_files.size() || hierarchy.size() != config.group_bin_files.size() )
        {
            std::cerr << "Hierarchy does not match with the number of provided files" << std::endl;
            return false;
        }
        else
        {
            for ( uint16_t h = 0; h < hierarchy.size(); ++h )
                config.filters[hierarchy[h]].push_back(
                    std::make_tuple( config.bloom_filter_files[h], config.group_bin_files[h] ) );
        }
    }
    return true;
}

} // namespace detail


bool GanonClassify::run( Config config )
{
    // disable on testing, messing up syncing with ctest
    if ( !config.testing )
        std::ios_base::sync_with_stdio( false ); // speed up output to STDOUT (allows buffering) ->
                                                 // https://en.cppreference.com/w/cpp/io/ios_base/sync_with_stdio

    // Check parameters
    config.clas_threads = config.threads - 2; //-1 reading, -1 printing
    if ( !config.output_unclassified_file.empty() )
    {
        config.output_unclassified = true;
        config.clas_threads        = config.clas_threads - 1; //-1 printing unclassified
    }
    else
    {
        config.output_unclassified = false;
    }

    if ( config.max_error_unique >= 0 && config.max_error_unique < config.max_error )
    {
        config.unique_filtering = true;
    }
    else
    {
        config.unique_filtering = false;
    }
    if ( !detail::parse_hierarchy( config ) )
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
    if ( !config.output_file.empty() )
    { // output to a file
        out.open( config.output_file );
    }
    else
    {
        out.copyfmt( std::cout ); // STDOUT
        out.clear( std::cout.rdstate() );
        out.basic_ios< char >::rdbuf( std::cout.rdbuf() );
    }

    // Queues for internal read handling
    // queue1 get reads from file
    // queue2 will get unclassified reads if hierachy == 2
    // if hierachy == 3 queue1 is used for unclassified and so on
    SafeQueue< detail::ReadBatches > queue1;
    SafeQueue< detail::ReadBatches > queue2;

    // Queues for classified, unclassified reads (print)
    SafeQueue< detail::ReadOut > classified_reads_queue;
    SafeQueue< detail::ReadOut > unclassified_reads_queue;

    // Control flags
    bool finished_reading     = false;
    bool finished_classifying = false;

    // Statistics values
    detail::Stats stats;

    // num_of_batches*num_of_reads_per_batch = max. amount of reads in memory
    int num_of_batches         = 1000;
    int num_of_reads_per_batch = 400;

    std::vector< std::future< void > > read_write;

    // Thread for reading input files
    timeLoadReads.start();
    read_write.emplace_back( std::async( std::launch::async, [=, &queue1, &finished_reading, &timeLoadReads, &stats] {
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
                // std::cerr << queue1->size() << std::endl;
                while ( queue1.size() > num_of_batches )
                {
                    ; // spin
                }
                seqan::StringSet< seqan::CharString > ids;
                seqan::StringSet< seqan::Dna5String > seqs;
                seqan::readRecords( ids, seqs, seqFileIn, num_of_reads_per_batch );
                stats.totalReads += seqan::length( ids );
                queue1.push( detail::ReadBatches{ ids, seqs } );
            }
            seqan::close( seqFileIn );
        }
        finished_reading = true;
        timeLoadReads.end();
    } ) );

    // Thread for printing classified reads
    timePrintClass.start();
    read_write.emplace_back(
        std::async( std::launch::async, [=, &classified_reads_queue, &out, &finished_classifying, &timePrintClass] {
            while ( true )
            {
                detail::ReadOut ro = classified_reads_queue.pop();
                for ( uint32_t i = 0; i < ro.matches.size(); ++i )
                    out << ro.readID << '\t' << ro.matches[i].group << '\t' << ro.matches[i].kmer_count << '\n';
                if ( finished_classifying && classified_reads_queue.empty() )
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
        read_write.emplace_back(
            std::async( std::launch::async,
                        [=, &unclassified_reads_queue, &out_unclassified, &finished_classifying, &timePrintUnclass] {
                            while ( true )
                            {
                                detail::ReadOut rou = unclassified_reads_queue.pop();
                                if ( rou.readID != "" ) // if not empty
                                    out_unclassified << rou.readID << '\n';
                                if ( finished_classifying && unclassified_reads_queue.empty() )
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
    uint16_t hierarchy_size = config.filters.size();
    // For every hiearchy level
    for ( auto const& hierarchy : config.filters )
    {
        ++hierarchy_id;
        // filter[h] = vector<(filter, map)> -> map already sort by key

        // std::string -> hierarchy.first [hierarchy_name]
        // std::vector< std::tuple< std::string, std::string > > = hierachy.second [hierarchy_files]

        timeLoadFilters.start();
        std::vector< detail::Filter > filter_hierarchy;
        detail::load_filters( filter_hierarchy, hierarchy.second );
        timeLoadFilters.end();

        // Exchange queues instance pointers for each hierachy (if not first)
        if ( hierarchy_id > 1 )
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

        // Threads for classification
        timeClass.start();
        for ( uint16_t taskNo = 0; taskNo < config.clas_threads; ++taskNo )
        {
            tasks.emplace_back( std::async(
                std::launch::async,
                [=,
                 &filter_hierarchy,
                 &classified_reads_queue,
                 &unclassified_reads_queue,
                 &finished_reading,
                 &stats,
                 &config] {
                    while ( true )
                    {
                        // std::cerr << pointer_current.size() << std::endl; //check if queue is getting empty (print 0's)
                        detail::ReadBatches rb = pointer_current->pop();
                        if ( rb.ids != "" )
                        {
                            detail::ReadBatches left_over_reads; // store unclassified reads for next iteration
                            for ( uint32_t readID = 0; readID < seqan::length( rb.ids ); ++readID )
                            {
                                uint16_t readLen = seqan::length( rb.seqs[readID] );
                                // count lens just once
                                if ( hierarchy_id == 1 )
                                    stats.sumReadLen += readLen;

                                // k-mer sizes should be the same among filters, groups should not overlap
                                uint16_t kmerSize  = filter_hierarchy[0].kmerSize;
                                uint16_t threshold = detail::kmer_threshold( readLen, kmerSize, config.max_error );

                                // Classify reads, returing matches and max kmer count achieved for the read
                                detail::Tmatches matches;
                                uint16_t         maxKmerCountRead =
                                    detail::classify_read( matches, filter_hierarchy, rb.seqs[readID], threshold );

                                // Filter reads by number of errors, special filtering for unique matches
                                detail::ReadOut read_out( rb.ids[readID] );
                                uint32_t        count_filtered_matches = detail::filter_matches(
                                    read_out, matches, kmerSize, readLen, maxKmerCountRead, config );

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

                        if ( finished_reading && pointer_current->empty() )
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
        timeClass.end();
    }
    finished_classifying = true;

    for ( auto&& task : read_write )
    {
        task.get();
    }

    if ( config.output_file != "" )
        out.close();

    if ( config.output_unclassified )
        out_unclassified.close();

    timeGanon.end();

    std::cerr << std::endl;
    detail::print_time(
        config, timeGanon, timeLoadReads, timeLoadFilters, timeClass, timePrintClass, timePrintUnclass );
    detail::print_stats( stats, timeClass );

    return true;
}
