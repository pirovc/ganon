# ===========================================================================
# ganon
# ===========================================================================

cmake_minimum_required( VERSION 3.2 FATAL_ERROR )
project( ganon LANGUAGES CXX )

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

set( CMAKE_CXX_STANDARD 17 )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
endif()

add_compile_options( -Wuninitialized -W -Wall -Wstrict-aliasing -pedantic
    -Wno-long-long -Wno-variadic-macros -Wunused -msse4.2 -O3 -DNDEBUG
    -static -march=native )

option( VERBOSE_CONFIG "Verbose mode for quick build setup debugging" OFF )

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# 1. threads:
find_package( Threads REQUIRED )

# 2. cxxopts:

add_library( cxxopts INTERFACE )
target_include_directories( cxxopts INTERFACE libs/cxxopts/include )

# 3. sdsl-lite:

add_library( sdsl-lite INTERFACE )
target_include_directories( sdsl-lite INTERFACE libs/sdsl-lite/include )

# 4. Zlib (optional for SeqAn):

find_package( ZLIB )

# 5. SeqAn:

set( SEQAN_INCLUDE_PATH libs/seqan/include )
find_package( SeqAn CONFIG REQUIRED PATHS libs/seqan/util/cmake )

add_library( seqan INTERFACE )
target_include_directories( seqan INTERFACE ${SEQAN_INCLUDE_DIRS} )
target_link_libraries( seqan INTERFACE ${SEQAN_LIBRARIES} )

add_compile_options( ${SEQAN_DEFINITIONS} )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}" )

# ----------------------------------------------------------------------------
# Verbose log
# ----------------------------------------------------------------------------

if( VERBOSE_CONFIG )
    message( STATUS "SeqAn symbols")
    message( STATUS "  SEQAN_DEFINITIONS : ${SEQAN_DEFINITIONS}" )
    message( STATUS "  SEQAN_CXX_FLAGS   : ${SEQAN_CXX_FLAGS}" )
    message( STATUS "  SEQAN_INCLUDE_DIRS: ${SEQAN_INCLUDE_DIRS}" )
    message( STATUS "  SEQAN_LIBRARIES   : ${SEQAN_LIBRARIES}" )
    message( STATUS "Misc symbols")
    message( STATUS "  CMAKE_CXX_FLAGS   : ${CMAKE_CXX_FLAGS}" )

    get_directory_property( dirCompileOptions COMPILE_OPTIONS )
    message( STATUS "  COMPILE_OPTIONS   : ${dirCompileOptions}" )
endif()

# ----------------------------------------------------------------------------
# Directories
# ----------------------------------------------------------------------------

add_subdirectory( src )
