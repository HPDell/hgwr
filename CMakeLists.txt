cmake_minimum_required(VERSION 3.12.0)
project(libhgwr VERSION 0.3.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/test/data)

add_subdirectory(src)

include(CTest)
enable_testing()

add_subdirectory(test)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
