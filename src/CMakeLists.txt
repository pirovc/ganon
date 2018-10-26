# ===========================================================================
# ganon utils
# ===========================================================================

# utils lib: this target holds resources shared by all other targets.

set( utils-src
    utils/include/utils/safequeue.hpp )

add_library( utils INTERFACE )
target_include_directories( utils INTERFACE utils/include )

# helper function:

function( link_to_ganon_dependencies target )
    target_link_libraries( ${target}
        PRIVATE
        sdsl-lite
        seqan
        Threads::Threads
        cxxopts
        utils )
endfunction()

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
