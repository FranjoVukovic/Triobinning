cmake_minimum_required(VERSION 3.15)
project(Triobinning VERSION 1.0)
add_executable(Triobinning main.cpp)

include(CTest)
set(CMAKE_CXX_STANDARD 17)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC external/bioparser/include)
target_include_directories(${PROJECT_NAME} PUBLIC external/thread_pool/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

#find_package(Git QUIET)
#if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
#    option(GIT_SUBMODULE "Check submodules during build" ON)
#    if(GIT_SUBMODULE)
#        message(STATUS "Submodule update")
#        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
#                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
#        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#        endif()
#    endif()
#endif()

#if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/bioparser/CMakeLists.txt")
#    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
#endif()

#if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/thread_pool/CMakeLists.txt")
#    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
#endif()

add_subdirectory(external/bioparser)
add_subdirectory(external/thread_pool)
add_subdirectory(external/googletest)
add_subdirectory(tests)

target_link_libraries(Triobinning PRIVATE z)