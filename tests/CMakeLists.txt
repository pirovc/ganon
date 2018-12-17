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

file(
    COPY
    ganon-classify/data/bacteria_group_bin.txt
    ganon-classify/data/classify.filter
    ganon-classify/data/classify_output.txt
    ganon-classify/data/simulated.1.fq
    ganon-classify/data/simulated.2.fq
    DESTINATION
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

add_executable( ganon-tests
    aux/Aux.hpp
    setup/Catch2.setup.cpp
    ganon-build/GanonBuild.test.cpp
    ganon-classify/GanonClassify.test.cpp
    utils/SafeQueue.test.cpp)

target_link_libraries( ganon-tests
    PRIVATE
    utils
    catch2
    ganon-build-lib
    ganon-classify-lib )

target_include_directories( ganon-tests PRIVATE . )

add_test( NAME GanonTestSuite
    COMMAND $<TARGET_FILE:ganon-tests>
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
