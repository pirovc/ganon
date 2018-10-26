#pragma once

#include <mutex>
#include <queue>

template < class T >
class SafeQueue
{
public:
    void push( T t )
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        m_queue.push( t );
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

    int size()
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        return m_queue.size();
    }

    bool empty()
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        return m_queue.empty();
    }

private:
    std::queue< T > m_queue;
    std::mutex      m_mutex;
};
