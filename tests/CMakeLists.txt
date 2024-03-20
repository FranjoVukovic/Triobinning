if (BUILD_TESTING)
    find_package(ZLIB REQUIRED)

    add_executable(tests my_test.cpp)
    target_link_libraries(tests PRIVATE GTest::gtest_main)
    target_link_libraries(tests PRIVATE ZLIB::ZLIB)

    include(GoogleTest)
    gtest_discover_tests(tests)
endif()