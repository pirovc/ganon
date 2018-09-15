# ===========================================================================
#                           ganon
# ===========================================================================

cmake_minimum_required( VERSION 3.0.0 )
project( ganon CXX )
message( STATUS "Configuring ganon" )

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package( ZLIB )
find_package( SeqAn CONFIG REQUIRED PATHS libs/seqan/util/cmake )

add_library( sdsl-lite INTERFACE IMPORTED )
target_include_directories( sdsl-lite INTERFACE libs/sdsl-lite/include )

add_library( seqan INTERFACE IMPORTED )
target_include_directories( seqan INTERFACE libs/seqan/include )

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

add_definitions( ${SEQAN_DEFINITIONS} )

set( CMAKE_CXX_FLAGS "-std=c++1z -Wuninitialized -W -Wall -Wstrict-aliasing -pedantic -Wno-long-long -Wno-variadic-macros -Wunused -msse4.2 -O3 -DNDEBUG -static -pthread -march=native" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}" )

# ----------------------------------------------------------------------------
# Build Targets
# ----------------------------------------------------------------------------

# Update the list of file names below if you add source files to your application.
add_executable( ganon-build ganon-build.cpp )
add_executable( ganon-update ganon-update.cpp )
add_executable( ganon-classify ganon-classify.cpp )

target_link_libraries( ganon-build sdsl-lite seqan )
target_link_libraries( ganon-update sdsl-lite seqan )
target_link_libraries( ganon-classify sdsl-lite seqan )
