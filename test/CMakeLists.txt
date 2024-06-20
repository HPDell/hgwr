cmake_minimum_required(VERSION 3.12.0)

include_directories(../include)

find_package(Catch2 3 REQUIRED)

if(NOT MSVC)
    add_definitions(-DDBL_MAX=__DBL_MAX__)
endif()

set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
add_definitions(-DTEST_DATA_DIR="${TEST_DATA_DIR}")

add_executable(testHGWRBML test_hgwrbfml.cpp helper.cpp)
target_link_libraries(testHGWRBML PRIVATE hgwrbml Catch2::Catch2WithMain)
add_test(
    NAME Test_HGWRBML
    COMMAND $<TARGET_FILE:testHGWRBML>  --success
)
