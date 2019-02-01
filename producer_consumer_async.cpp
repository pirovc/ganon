/**
 * g++ -std=c++11 -pthread producer_consumer_asyn.cpp
 **/

#include <iostream>                     
#include <mutex>              
#include <condition_variable> 
#include <unistd.h>
#include <future>
#include "src/utils/include/utils/SafeQueue.hpp"



int main(){

  std::condition_variable cv;
  std::mutex mtx;

  SafeQueue< std::string > queue;
  int total_reads = 100;
  bool finished = false;
  int threads = 4;

  std::future< void > read_task( 
    std::async( std::launch::async, [=, &queue, &finished , &cv] 
    {
        for( int r = 0; r<total_reads;++r){
          queue.push(std::to_string(r));
          usleep(50000);
          //std::cout << r << " added to the queue" << std::endl;
          cv.notify_one(); // every read added, notify one waiting thread
        }
        finished = true;
        cv.notify_all(); // notify all to exit
    } 
    ) 
  );


  std::vector< std::future< void > > tasks;
  for ( int taskNo = 0; taskNo < threads; ++taskNo )
  {
      tasks.emplace_back( 
        std::async( std::launch::async, [=, &queue, &finished, &cv, &mtx] 
        {
            
            bool over = false;
            while(true){

              while ( queue.empty() ){
                if(finished){
                  over = true;
                  break;
                }
                std::cout << "Waiting on thread "  << taskNo << std::endl;
                std::unique_lock<std::mutex> lck(mtx);
                cv.wait(lck);
              }
              if (over) break;

              std::cout << queue.size() << std::endl;
              std::string r = queue.pop();
              if ( r != "" )
              {
                  std::cout << "Read " << r << " on the thread "  << taskNo << std::endl;
              }

              
            }
            
            
            
        }
        )
      );
  }

  read_task.get();
  for ( auto&& task : tasks )
  {
      task.get();
  }

  return 0;
}
