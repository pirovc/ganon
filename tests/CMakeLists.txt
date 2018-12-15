# ===========================================================================
# ganon tests
# ===========================================================================

file(
    COPY
        ganon-build/data/bacteria_NC_010333.1.fasta.gz
        ganon-build/data/bacteria_NC_017163.1.fasta.gz
        ganon-build/data/bacteria_NC_017164.1.fasta.gz
        ganon-build/data/bacteria_NC_017543.1.fasta.gz
        ganon-build/data/bacteria_seqid_bin.txt
        ganon-build/data/build_output.filter
    DESTINATION
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

add_executable( ganon-tests
    setup/Catch2.setup.cpp
    ganon-build/GanonBuild.test.cpp
    utils/SafeQueue.test.cpp )

target_link_libraries( ganon-tests
    PRIVATE
    utils
    catch2
    ganon-build-lib )

# TODO resolve relative path to test target
# add_test( NAME GanonTestSuite COMMAND ganon-tests )
