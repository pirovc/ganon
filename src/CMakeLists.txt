# ===========================================================================
# ganon utils
# ===========================================================================

# utils lib: this target holds resources shared by all other targets.

set( utils-src
    utils/include/utils/safequeue.hpp )

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

set( ganon-build-src
    ganon-build/ganon-build.cpp )

add_executable( ganon-build ${ganon-build-src} )
link_to_ganon_dependencies( ganon-build )

# ----------------------------------------------------------------------------
# ganon-classify
# ----------------------------------------------------------------------------

set( ganon-classify-src
    ganon-classify/ganon-classify.cpp )

add_executable( ganon-classify ${ganon-classify-src} )
link_to_ganon_dependencies( ganon-classify )
