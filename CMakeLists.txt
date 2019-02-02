# =============================================================================
# ganon
# =============================================================================

cmake_minimum_required( VERSION 3.8.2 FATAL_ERROR )
project( ganon VERSION 0.9.0 LANGUAGES CXX )

# -----------------------------------------------------------------------------
# build setup
# -----------------------------------------------------------------------------

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

if( NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" )
    message( FATAL_ERROR
        "Compiler id '${CMAKE_CXX_COMPILER_ID}' is not supported, please \
        check the documentation." )
endif()

option( VERBOSE_CONFIG "Verbose mode for quick build setup debugging" OFF )

# -----------------------------------------------------------------------------
# build types
# -----------------------------------------------------------------------------

get_property( isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG )
if( isMultiConfig )
    if( NOT "Coverage" IN_LIST CMAKE_CONFIGURATION_TYPES )
        list( APPEND CMAKE_CONFIGURATION_TYPES Coverage )
    endif()
else()
    set( allowableBuildTypes Debug Release RelWithDebInfo MinSizeRel Coverage )

    set_property( CACHE CMAKE_BUILD_TYPE
        PROPERTY STRINGS "${allowableBuildTypes}" )

    if( NOT CMAKE_BUILD_TYPE )
        set( CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE )
    elseif( NOT CMAKE_BUILD_TYPE IN_LIST allowableBuildTypes )
        message( FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}" )
    endif()
endif()

# Release flags:

set( CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" )

# Coverage flags:

set( CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} -O0 --coverage" )

set( CMAKE_EXE_LINKER_FLAGS_COVERAGE
    "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -lgcov --coverage" )

set( CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" )
set( CMAKE_STATIC_LINKER_FLAGS_COVERAGE "${CMAKE_STATIC_LINKER_FLAGS_DEBUG}" )
set( CMAKE_MODULE_LINKER_FLAGS_COVERAGE "${CMAKE_MODULE_LINKER_FLAGS_DEBUG}" )

# warning flags:

add_compile_options( -Wall -Wextra -Wshadow -Wuninitialized -Wnon-virtual-dtor
    -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic
    -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion
    -Wformat=2 -Wstrict-aliasing -Wno-long-long -Wno-variadic-macros )

if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    add_compile_options( -Wmisleading-indentation -Wduplicated-cond
        -Wduplicated-branches -Wlogical-op -Wuseless-cast )
endif()

# seqan specific flags:

#add_compile_options( -static -msse4.2 -march=native )
add_compile_options( -static -march=native )

# -----------------------------------------------------------------------------
# dependencies and 3rd party libraries
# -----------------------------------------------------------------------------

# 1. threads:

find_package( Threads REQUIRED )

# 2. cxxopts:

add_library( cxxopts INTERFACE )
target_include_directories( cxxopts SYSTEM INTERFACE libs/cxxopts/include )

# 3. sdsl-lite:

add_library( sdsl-lite INTERFACE )
target_include_directories( sdsl-lite
    SYSTEM INTERFACE libs/sdsl-lite/include )

# 4. Zlib (optional for SeqAn):

find_package( ZLIB )

# 5. SeqAn:

set( SEQAN_INCLUDE_PATH libs/seqan/include )
find_package( SeqAn CONFIG REQUIRED PATHS libs/seqan/util/cmake )

add_library( seqan INTERFACE )
target_include_directories( seqan SYSTEM INTERFACE ${SEQAN_INCLUDE_DIRS} )
target_link_libraries( seqan INTERFACE ${SEQAN_LIBRARIES} )

add_compile_options( ${SEQAN_DEFINITIONS} )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}" )

# 6. Catch2:

add_library( catch2 INTERFACE )
target_include_directories( catch2
    SYSTEM INTERFACE libs/Catch2/single_include )

# -----------------------------------------------------------------------------
# verbose log
# -----------------------------------------------------------------------------

if( VERBOSE_CONFIG )
    message( STATUS "SeqAn symbols")
    message( STATUS "  SEQAN_DEFINITIONS : ${SEQAN_DEFINITIONS}" )
    message( STATUS "  SEQAN_CXX_FLAGS   : ${SEQAN_CXX_FLAGS}" )
    message( STATUS "  SEQAN_INCLUDE_DIRS: ${SEQAN_INCLUDE_DIRS}" )
    message( STATUS "  SEQAN_LIBRARIES   : ${SEQAN_LIBRARIES}" )
    message( STATUS "Misc symbols")
    message( STATUS "  Build type        : ${CMAKE_BUILD_TYPE}" )
    message( STATUS "  CMAKE_CXX_FLAGS   : ${CMAKE_CXX_FLAGS}" )

    get_directory_property( dirCompileOptions COMPILE_OPTIONS )
    message( STATUS "  COMPILE_OPTIONS   : ${dirCompileOptions}" )
endif()

# -----------------------------------------------------------------------------
# folders
# -----------------------------------------------------------------------------

enable_testing()
add_subdirectory( src )
add_subdirectory( tests )
