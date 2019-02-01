/**
 * Producer-Comsumer example, written in C++ May 4, 2014
 * Compiled on OSX 10.9, using:
 * g++ -std=c++11 producer_consumer.cpp
 * https://austingwalters.com/multithreading-producer-consumer-problem/
 * https://austingwalters.com/multithreading-semaphores/
 * https://austingwalters.com/mutex-process-synchronization/
 * https://www.justsoftwaresolutions.co.uk/threading/locks-mutexes-semaphores.html
 **/

#include <iostream>           
#include <thread>             
#include <mutex>              
#include <condition_variable> 

std::mutex mtx;
std::condition_variable cv;

int meal = 0;

/* Consumer */
void waiter(int ordernumber){
  std::unique_lock<std::mutex> lck(mtx);
  while(meal == 0) cv.wait(lck);
  std::cout << "Order: ";
  std::cout << ordernumber + 1 << " being taken care of with ";
  std::cout << meal - 1 << " meals also ready." << std::endl;
  meal--;
}

/* Producer */
void makeMeal(int ordernumber){
  std::unique_lock<std::mutex> lck(mtx);
  meal++;
  cv.notify_one();
}

int main(){

  std::thread chefs[10];
  std::thread waiters[10];

  /* Initialize customers and cheifs */
  for (int order = 0; order < 10; order++){
    chefs[order] = std::thread(makeMeal, order);
    waiters[order] = std::thread(waiter, order);
  }

  /* Join the threads to the main threads */
  for (int order = 0; order < 10; order++) {
    waiters[order].join();   
    chefs[order].join(); 
  }

  return 0;
}