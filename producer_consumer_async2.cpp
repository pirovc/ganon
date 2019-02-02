/**
 * g++ -std=c++11 -pthread producer_consumer_asyn.cpp
 **/

#include "src/utils/include/utils/SafeQueue.hpp"
#include <future>
#include <iostream>
#include <mutex>
#include <unistd.h>

std::mutex m;
int        main()
{


    SafeQueue< std::string > queue( 10 );
    int                      total_reads = 1000;
    bool                     finished    = false;
    int                      threads     = 48;

    std::future< void > read_task( std::async( std::launch::async, [=, &queue] {
        for ( int r = 0; r < total_reads; ++r )
        {
            queue.push( std::to_string( r ) );
            std::cout << "Read " << r << " pushed" << std::endl;
            // usleep(50000);
        }
        queue.notify_push_over();
    } ) );


    std::vector< std::future< void > > tasks;
    for ( int taskNo = 0; taskNo < threads; ++taskNo )
    {
        tasks.emplace_back( std::async( std::launch::async, [=, &queue] {

            while ( true )
            {

                std::string r = queue.pop();
                if ( r != "" )
                {
                    std::cout << "Read " << r << " poped on the thread " << taskNo << std::endl;
                    usleep( 50000 );
                }
                else
                {
                    break;
                }
            }


        } ) );
    }

    read_task.get();
    for ( auto&& task : tasks )
    {
        task.get();
    }

    return 0;
}
