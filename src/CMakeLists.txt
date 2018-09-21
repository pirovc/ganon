# utils:

set( utils-src
    utils/safequeue.hpp )

add_library( utils ${utils-src} )
set_target_properties( utils PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 17 )
target_include_directories( utils PUBLIC utils )

# ganon-build:

set( ganon-build-src
    ganon-build/ganon-build.cpp )

add_executable( ganon-build ${ganon-build-src} )

# ganon-classify:

set( ganon-classify-src
    ganon-classify/ganon-classify.cpp )

add_executable( ganon-classify ${ganon-classify-src} )

# ganon-update:

set( ganon-update-src
    ganon-update/ganon-update.cpp )

add_executable( ganon-update ${ganon-update-src} )

set( allTargets ganon-build ganon-classify ganon-update )

foreach( target ${allTargets} )
    target_link_libraries(
        ${target}
        PRIVATE
        sdsl-lite
        seqan
        Threads::Threads
        cxxopts
        utils )
endforeach()
