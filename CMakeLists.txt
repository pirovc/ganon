# ===========================================================================
#                           ganon
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (ganon CXX)
message (STATUS "Configuring ganon")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package (ZLIB)
find_package (SeqAn CONFIG REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SDSL_INCLUDE_DIRS})
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable( ganon-build ganon-build.cpp)
add_executable( ganon-update ganon-update.cpp)
add_executable( ganon-classify ganon-classify.cpp)

target_link_libraries (ganon-build ${SEQAN_LIBRARIES})
target_link_libraries (ganon-update ${SEQAN_LIBRARIES})
target_link_libraries (ganon-classify ${SEQAN_LIBRARIES})
