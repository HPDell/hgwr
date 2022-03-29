cmake_minimum_required(VERSION 3.12.0)

include_directories(../include)

find_package(Armadillo REQUIRED)
if(ARMADILLO_FOUND)
    add_definitions(-DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_DONT_USE_WRAPPER)
    include_directories(${ARMADILLO_INCLUDE_DIRS})

    find_package(BLAS)
    if(NOT BLAS_FOUND)
        message(SEND_ERROR "BLAS Not Found!")
    endif(NOT BLAS_FOUND)
    
    find_package(LAPACK)
    if(NOT LAPACK_FOUND)
        message(SEND_ERROR "LAPACK Not Found!")
    endif(NOT LAPACK_FOUND)
else(ARMADILLO_FOUND)
    message(FATAL_ERROR "Armadillo Not Found!")
endif(ARMADILLO_FOUND)

find_package(Boost COMPONENTS program_options)
if(Boost_FOUND)
    message(STATUS "Boost found")
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost Library Cannot Found!")
endif()

set(ADDON_DEBUG_LIBRARIES
    ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG} 
)

set(ADDON_OPTIMIZED_LIBRARIES
    ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
)

add_executable(testHGWRBML main.cpp)
target_link_libraries(testHGWRBML PRIVATE hgwrbml debug ${ADDON_DEBUG_LIBRARIES} optimized ${ADDON_OPTIMIZED_LIBRARIES})

add_test(
    NAME Test_HGWRBML
    COMMAND $<TARGET_FILE:testHGWRBML> -d ${TEST_DATA_DIR} -b 64
    WORKING_DIRECTORY $<TARGET_FILE_DIR:testHGWRBML>
)