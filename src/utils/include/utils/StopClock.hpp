#pragma once

#include <chrono>
#include <ctime>
#include <iomanip>

class StopClock
{
public:
    using Clock     = std::chrono::system_clock;
    using TimePoint = std::chrono::time_point< Clock >;
    using Seconds   = double;

    void start()
    {
        m_beginRound = Clock::now();

        if ( m_firstStart )
        {
            m_begin      = m_beginRound;
            m_firstStart = false;
        }
    }

    void stop()
    {
        m_end = Clock::now();
        m_runTime += m_end - m_beginRound;
    }

    Seconds elapsed() const noexcept
    {
        return m_runTime.count();
    }

    TimePoint begin() const noexcept
    {
        return m_begin;
    }

    TimePoint end() const noexcept
    {
        return m_end;
    }

private:
    bool                             m_firstStart{ true };
    TimePoint                        m_begin;
    TimePoint                        m_beginRound;
    TimePoint                        m_end;
    std::chrono::duration< Seconds > m_runTime{ 0.0 };
};

template < typename Stream >
inline Stream& operator<<( Stream& stream, const StopClock::TimePoint& timepoint )
{
    const auto time = std::chrono::system_clock::to_time_t( timepoint );
    stream << std::put_time( std::localtime( &time ), "%c" );

    return stream;
}
