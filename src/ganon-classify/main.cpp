#include "CommandLineParser.hpp"
#include "GanonClassify.hpp"

#include <cstdlib>
#include <utility>

int main( int argc, char** argv )
{
    if ( auto config = CommandLineParser::parse( argc, argv ); config.has_value() )
    {
        return GanonClassify::run( std::move( config.value() ) ) ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
