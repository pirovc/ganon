# ===========================================================================
# ganon utils
# ===========================================================================

# utils lib: this target holds resources shared by all other targets.

set( utils-src
    utils/include/utils/SafeQueue.hpp
    utils/include/utils/Time.hpp )

add_library( utils STATIC ${utils-src} )
set_target_properties( utils PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 17 )
target_include_directories( utils PUBLIC utils/include )

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

set( ganon-classify-src
    ganon-classify/CommandLineParser.cpp
    ganon-classify/CommandLineParser.hpp
    ganon-classify/Config.hpp
    ganon-classify/GanonClassify.cpp
    ganon-classify/GanonClassify.hpp
    ganon-classify/main.cpp )

add_executable( ganon-classify ${ganon-classify-src} )
link_to_ganon_dependencies( ganon-classify )
