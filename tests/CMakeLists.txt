# ===========================================================================
# ganon tests
# ===========================================================================

add_executable( ganon-tests
    setup/Catch2.setup.cpp
    utils/SafeQueue.test.cpp )

target_link_libraries( ganon-tests
    PRIVATE
    utils
    catch2 )

add_test( NAME GanonTestSuite COMMAND ganon-tests )
