cmake_minimum_required(VERSION 3.12.0)
project(libhgwr VERSION 0.3.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

option(WITH_TESTS "Whether to build tests" ON)

add_subdirectory(src)

if(WITH_TESTS)
    include(CTest)
    set(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/test/data)
    enable_testing()
    add_subdirectory(test)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
