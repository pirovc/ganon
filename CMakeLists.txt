# ===========================================================================
#                           ganon
# ===========================================================================

cmake_minimum_required( VERSION 3.2 FATAL_ERROR )
project( ganon CXX )

option( VERBOSE_CONFIG "Verbose mode for build setup" OFF )

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# 1. sdsl-lite:

add_library( sdsl-lite INTERFACE )
target_include_directories( sdsl-lite INTERFACE libs/sdsl-lite/include )

# 2. SeqAn:

find_package( ZLIB ) # optional for SeqAn

set( SEQAN_INCLUDE_PATH libs/seqan/include )
find_package( SeqAn CONFIG REQUIRED PATHS libs/seqan/util/cmake )

add_library( seqan INTERFACE )
target_include_directories( seqan INTERFACE ${SEQAN_INCLUDE_DIRS} )
target_link_libraries( seqan INTERFACE ${SEQAN_LIBRARIES} )

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

add_definitions( ${SEQAN_DEFINITIONS} )

set( CMAKE_CXX_FLAGS "-std=c++1z -Wuninitialized -W -Wall -Wstrict-aliasing -pedantic -Wno-long-long -Wno-variadic-macros -Wunused -msse4.2 -O3 -DNDEBUG -static -pthread -march=native" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}" )

if( VERBOSE_CONFIG )
    message( STATUS "SeqAn symbols")
    message( STATUS "  SEQAN_DEFINITIONS : ${SEQAN_DEFINITIONS}" )
    message( STATUS "  SEQAN_CXX_FLAGS   : ${SEQAN_CXX_FLAGS}" )
    message( STATUS "  SEQAN_INCLUDE_DIRS: ${SEQAN_INCLUDE_DIRS}" )
    message( STATUS "  SEQAN_LIBRARIES   : ${SEQAN_LIBRARIES}" )
    message( STATUS "Misc symbols")
    message( STATUS "  CMAKE_CXX_FLAGS   : ${CMAKE_CXX_FLAGS}" )
endif()

# ----------------------------------------------------------------------------
# Build Targets
# ----------------------------------------------------------------------------

add_executable( ganon-build ganon-build.cpp )
add_executable( ganon-update ganon-update.cpp )
add_executable( ganon-classify ganon-classify.cpp )

target_link_libraries( ganon-build sdsl-lite seqan )
target_link_libraries( ganon-update sdsl-lite seqan )
target_link_libraries( ganon-classify sdsl-lite seqan )
