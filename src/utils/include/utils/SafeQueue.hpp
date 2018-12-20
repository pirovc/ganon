#pragma once

#include <mutex>
#include <optional>
#include <queue>
#include <type_traits>
#include <utility>

template < class T,
           class Enabled = std::enable_if_t< std::is_move_assignable< T >::value && //
                                             std::is_move_constructible< T >::value > >
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

    template < class... Args >
    void emplace( Args&&... args )
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        m_queue.emplace( std::forward< Args >( args )... );
    }

    std::optional< T > pop()
    {
        std::lock_guard< std::mutex > lock( m_mutex );
        if ( m_queue.empty() )
        {
            return std::nullopt;
        }

        auto value = std::make_optional( std::move( m_queue.front() ) );
        m_queue.pop();
        return value;
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
