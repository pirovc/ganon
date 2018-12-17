# ===========================================================================
# ganon utils
# ===========================================================================

add_library( utils INTERFACE )
target_include_directories( utils INTERFACE utils/include )

# ===========================================================================
# config files
# ===========================================================================

configure_file(
    ${PROJECT_SOURCE_DIR}/src/defaults/defaults.hpp.in
    ${PROJECT_BINARY_DIR}/include/defaults/defaults.hpp )

add_library( defaults INTERFACE )

target_include_directories( defaults 
    INTERFACE 
    ${PROJECT_BINARY_DIR}/include )

# ===========================================================================
# ganon targets
# ===========================================================================

# helper function:

function( link_to_ganon_dependencies target )
    target_link_libraries( ${target}
        PRIVATE
        sdsl-lite
        seqan
        Threads::Threads
        cxxopts
        utils
        defaults )
endfunction()

# ----------------------------------------------------------------------------
# ganon-build
# ----------------------------------------------------------------------------

set( ganon-build-lib-src
    ganon-build/CommandLineParser.cpp
    ganon-build/GanonBuild.cpp )

add_library( ganon-build-lib STATIC ${ganon-build-lib-src} )
link_to_ganon_dependencies( ganon-build-lib )

target_include_directories( ganon-build-lib
    PUBLIC
    ganon-build/include
    PRIVATE
    ganon-build/include/ganon-build )

add_executable( ganon-build ganon-build/main.cpp )
target_link_libraries( ganon-build ganon-build-lib )

# ----------------------------------------------------------------------------
# ganon-classify
# ----------------------------------------------------------------------------

set( ganon-classify-lib-src
    ganon-classify/CommandLineParser.cpp
    ganon-classify/GanonClassify.cpp )

add_library( ganon-classify-lib STATIC ${ganon-classify-lib-src} )
link_to_ganon_dependencies( ganon-classify-lib )

target_include_directories( ganon-classify-lib
    PUBLIC
    ganon-classify/include
    PRIVATE
    ganon-classify/include/ganon-classify )

add_executable( ganon-classify ganon-classify/main.cpp )
target_link_libraries( ganon-classify ganon-classify-lib )
