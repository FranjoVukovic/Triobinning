cmake_minimum_required(VERSION 3.15)

project(Triobinning VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FetchContent)

find_package(biosoup 0.11.0 QUIET)
if (NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.11.0)

  FetchContent_GetProperties(biosoup)
  if (NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(
      ${biosoup_SOURCE_DIR}
      ${biosoup_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(bioparser 3.0.15 QUIET)
if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.15)

  FetchContent_GetProperties(bioparser)
  if (NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(
      ${bioparser_SOURCE_DIR}
      ${bioparser_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(thread_pool 4.0.0 QUIET)
if (NOT thread_pool_FOUND)
  FetchContent_Declare(
    thread_pool
    GIT_REPOSITORY https://github.com/rvaser/thread_pool
    GIT_TAG 4.0.0)

  FetchContent_GetProperties(thread_pool)
  if (NOT thread_pool_POPULATED)
    FetchContent_Populate(thread_pool)
    add_subdirectory(
      ${thread_pool_SOURCE_DIR}
      ${thread_pool_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()


find_package(GTest 1.14.0 QUIET)
  if (NOT GTest_FOUND)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG v1.14.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()

add_library(Triobinning
  src/kmer_maker.cpp
  src/kmer_finder.cpp)

target_link_libraries(Triobinning
biosoup::biosoup
bioparser::bioparser
thread_pool::thread_pool)

target_include_directories(Triobinning PUBLIC
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

add_executable(triobinning_exe src/main.cpp)
target_link_libraries(triobinning_exe Triobinning)

target_compile_definitions(Triobinning PUBLIC TRIOBINNING_VERSION="${PROJECT_VERSION}")
set_property(TARGET triobinning_exe PROPERTY OUTPUT_NAME Triobinning)

install(TARGETS triobinning_exe DESTINATION ${CMAKE_INSTALL_BINDIR})

enable_testing()
add_executable(triobinning_test tests/my_test.cpp)
target_link_libraries(triobinning_test 
Triobinning 
bioparser::bioparser
GTest::Main)

gtest_discover_tests(triobinning_test)