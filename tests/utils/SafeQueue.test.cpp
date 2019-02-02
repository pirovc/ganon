#include <utils/SafeQueue.hpp>

#include <catch2/catch.hpp>

SCENARIO( "Pushing an element into an empty queue", "[utils][safequeue]" )
{
    GIVEN( "An empty queue" )
    {
        SafeQueue< int > queue;

        REQUIRE( queue.size() == 0 );
        REQUIRE( queue.empty() );

        WHEN( "one element is pushed" )
        {
            queue.push( 13 );
            queue.notify_push_over();

            THEN( "the queue size is one and the queue is not empty anymore" )
            {
                REQUIRE( queue.size() == 1 );
                REQUIRE_FALSE( queue.empty() );
            }
        }
    }
}

SCENARIO( "Popping from an empty queue", "[utils][safequeue]" )
{
    GIVEN( "An empty queue" )
    {
        SafeQueue< int > queue;
        queue.notify_push_over();

        REQUIRE( queue.size() == 0 );
        REQUIRE( queue.empty() );

        WHEN( "one element is popped" )
        {
            const auto value = queue.pop();

            THEN( "the popped element is the default value of the type" )
            {
                REQUIRE( value == int{} );
            }
            AND_THEN( "the queue remains empty" )
            {
                REQUIRE( queue.size() == 0 );
                REQUIRE( queue.empty() );
            }
        }
    }
}

SCENARIO( "FIFO property", "[utils][safequeue]" )
{
    GIVEN( "A queue with three elements" )
    {
        SafeQueue< int > queue;

        constexpr int firstPushedValue  = 11;
        constexpr int secondPushedValue = 12;
        constexpr int thirdPushedValue  = 13;

        queue.push( firstPushedValue );
        queue.push( secondPushedValue );
        queue.push( thirdPushedValue );
        queue.notify_push_over();

        REQUIRE( queue.size() == 3 );
        REQUIRE_FALSE( queue.empty() );

        WHEN( "the three elements are popped" )
        {
            const auto firstPoppedValue  = queue.pop();
            const auto secondPoppedValue = queue.pop();
            const auto thirdPoppedValue  = queue.pop();

            THEN( "the elements are popped on the same order they were pushed" )
            {
                REQUIRE( firstPoppedValue == firstPushedValue );
                REQUIRE( secondPoppedValue == secondPushedValue );
                REQUIRE( thirdPoppedValue == thirdPushedValue );
            }
            AND_THEN( "the queue is empty" )
            {
                REQUIRE( queue.size() == 0 );
                REQUIRE( queue.empty() );
            }
        }
    }
}
