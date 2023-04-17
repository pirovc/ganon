# =============================================================================
# ganon
# =============================================================================

cmake_minimum_required( VERSION 3.10 FATAL_ERROR )
project( ganon VERSION 1.6.0 LANGUAGES CXX )

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
option( CONDA "Flag for compilation in conda env." OFF )
option( LONGREADS "Uses uint32_t for count in ganon-classify. Useful for very long reads (>65535bp)" OFF )
option( INCLUDE_DIRS "Include directories to look for libraries" "" )

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

add_compile_options( -Wall -Wextra -Wshadow -Wuninitialized -Wcast-align -Wunused 
    -Woverloaded-virtual -Wpedantic -Wnull-dereference -Wdouble-promotion 
    -Wformat=2 -Wstrict-aliasing -Wunused-variable )

add_compile_options( -Wno-shadow -Wno-old-style-cast )

if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    add_compile_options( -Wmisleading-indentation -Wduplicated-cond
        -Wduplicated-branches -Wlogical-op -Wuseless-cast )
endif()

# seqan specific flags:

if ( NOT CONDA )
	add_compile_options( -static -march=native )
endif()

if( LONGREADS )
    add_compile_options(-DLONGREADS) 
endif()

# -----------------------------------------------------------------------------
# dependencies and 3rd party libraries
# -----------------------------------------------------------------------------

# 1. threads:

find_package( Threads REQUIRED )

# 2. cxxopts:

add_library( cxxopts INTERFACE )
if( INCLUDE_DIRS )
    target_include_directories( cxxopts SYSTEM INTERFACE ${INCLUDE_DIRS} )
else()
    target_include_directories( cxxopts SYSTEM INTERFACE libs/cxxopts/include )
endif()

# 3. Zlib (optional for SeqAn):

find_package( ZLIB )

# 4. SeqAn3:

set( SEQAN3_INCLUDE_PATH libs/seqan3/include )
find_package (SeqAn3 3.1.0 REQUIRED HINTS libs/seqan3/build_system)

add_library( seqan3 INTERFACE )
target_include_directories( seqan3 SYSTEM INTERFACE ${SEQAN3_INCLUDE_DIRS} )
target_link_libraries( seqan3 INTERFACE ${SEQAN3_LIBRARIES} )

add_compile_options( ${SEQAN3_DEFINITIONS} )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS}" )

# 5. Catch2:

add_library( Catch2 INTERFACE )
if( INCLUDE_DIRS )
    target_include_directories( Catch2 SYSTEM INTERFACE ${INCLUDE_DIRS} )
else()
    target_include_directories( Catch2 SYSTEM INTERFACE libs/Catch2/single_include )
endif()

# 6. robin-hood-hashing:

add_library( robin-hood-hashing INTERFACE )
if( INCLUDE_DIRS )
    target_include_directories( robin-hood-hashing SYSTEM INTERFACE ${INCLUDE_DIRS} )
else()
    target_include_directories( robin-hood-hashing SYSTEM INTERFACE libs/robin-hood-hashing/src/include )
endif()

# -----------------------------------------------------------------------------
# verbose log
# -----------------------------------------------------------------------------

if( VERBOSE_CONFIG )
    message( STATUS "SeqAn3 symbols")
    message( STATUS "  SEQAN3_VERSION   : ${SEQAN3_VERSION}" )
    message( STATUS "  SEQAN3_DEFINITIONS   : ${SEQAN3_DEFINITIONS}" )
    message( STATUS "  SEQAN3_CXX_FLAGS     : ${SEQAN3_CXX_FLAGS}" )
    message( STATUS "  SEQAN3_INCLUDE_DIRS  : ${SEQAN3_INCLUDE_DIRS}" )
    message( STATUS "  SEQAN3_LIBRARIES     : ${SEQAN3_LIBRARIES}" )
    message( STATUS "Misc symbols")
    message( STATUS "  Build type          : ${CMAKE_BUILD_TYPE}" )
    message( STATUS "  CMAKE_CXX_FLAGS     : ${CMAKE_CXX_FLAGS}" )
    message( STATUS "  INCLUDE_DIRS        : ${INCLUDE_DIRS}" )
    message( STATUS "  CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}" )
    message( STATUS "  CONDA               : ${CONDA}" )
    message( STATUS "  LONGREADS           : ${LONGREADS}" )
    get_directory_property( dirCompileOptions COMPILE_OPTIONS )
    message( STATUS "  COMPILE_OPTIONS   : ${dirCompileOptions}" )

endif()

# -----------------------------------------------------------------------------
# folders
# -----------------------------------------------------------------------------

enable_testing()
add_subdirectory( src )
add_subdirectory( tests )
