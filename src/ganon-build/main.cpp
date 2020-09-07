#include <ganon-build/CommandLineParser.hpp>
#include <ganon-build/GanonBuild.hpp>

#include <cstdlib>
#include <utility>

int main( int argc, char** argv )
{
    if ( auto config = GanonBuild::CommandLineParser::parse( argc, argv ); config.has_value() )
    {
        return GanonBuild::run( std::move( config.value() ) ) ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    else
    {
        return argc == 1 ? EXIT_FAILURE : EXIT_SUCCESS;
    }
}
