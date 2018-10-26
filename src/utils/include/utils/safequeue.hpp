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
