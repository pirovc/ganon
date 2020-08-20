#include <ganon-classify/CommandLineParser.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <cstdlib>
#include <utility>

int main( int argc, char** argv )
{
    if ( auto config = GanonClassify::CommandLineParser::parse( argc, argv ); config.has_value() )
    {
        return GanonClassify::run( std::move( config.value() ) ) ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    else
    {
        return argc == 1 ? EXIT_FAILURE : EXIT_SUCCESS;
    }
}
