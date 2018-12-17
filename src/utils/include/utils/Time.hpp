#pragma once

#include <chrono>
#include <ctime>
#include <string>
#include <vector>

struct Time
{

    typedef std::chrono::time_point< std::chrono::high_resolution_clock > Tpoint;
    typedef std::chrono::duration< double >                               Telapsed;
    // array of start and end times to account for multiple iterations
    std::vector< Tpoint > _start;
    std::vector< Tpoint > _end;

    void start()
    {
        _start.emplace_back( std::chrono::high_resolution_clock::now() );
    }

    void end()
    {
        _end.emplace_back( std::chrono::high_resolution_clock::now() );
    }

    double get_elapsed()
    {
        Telapsed elapsed{ 0.0 };
        for ( auto i = 0u; i < _start.size(); ++i )
        {
            elapsed += _end[i] - _start[i];
        }
        return elapsed.count();
    }

    std::string get_start_ctime()
    {
        auto start_ctime = std::chrono::system_clock::to_time_t( _start.front() );
        return std::ctime( &start_ctime );
    }
    std::string get_end_ctime()
    {
        auto end_ctime = std::chrono::system_clock::to_time_t( _end.back() );
        return std::ctime( &end_ctime );
    }
};
