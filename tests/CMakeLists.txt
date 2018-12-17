# ===========================================================================
# ganon tests
# ===========================================================================

file( COPY ganon-build/data DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
file( COPY ganon-classify/data DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

add_executable( ganon-tests
    aux/Aux.hpp
    setup/Catch2.setup.cpp
    ganon-build/GanonBuild.test.cpp
    ganon-classify/GanonClassify.test.cpp
    utils/SafeQueue.test.cpp )

target_link_libraries( ganon-tests
    PRIVATE
    utils
    catch2
    ganon-build-lib
    ganon-classify-lib )

target_include_directories( ganon-tests PRIVATE . )

add_test( NAME GanonTestSuite
    COMMAND $<TARGET_FILE:ganon-tests>
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data )
