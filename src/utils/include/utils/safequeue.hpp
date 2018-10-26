#pragma once

#include <mutex>
#include <queue>
#include <utility>

template < class T >
class SafeQueue
{
public:
    void push( const T& value )
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        m_queue.push( value );
    }

    void push( T&& value )
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        m_queue.push( value );
    }

    T pop()
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        if ( m_queue.empty() )
            return T();
        T val = m_queue.front();
        m_queue.pop();
        return val;
    }

    auto size() const
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        return m_queue.size();
    }

    bool empty() const
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        return m_queue.empty();
    }

private:
    std::queue< T >    m_queue;
    mutable std::mutex m_mutex;
};
