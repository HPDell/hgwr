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

find_package(OpenMP)
if(OpenMP_FOUND AND OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
    SET (ENABLE_OpenMP TRUE CACHE BOOL "Determines whether OpemMP support should be built")
    if(ENABLE_OpenMP)
        add_definitions(-DENABLE_OPENMP)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        include_directories(${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS})
    endif(ENABLE_OpenMP)
endif()

find_package(GSL)
if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "GSL cannot found!")
endif()

if(OpenMP_FOUND AND ENABLE_OpenMP)
    set(ADDON_LIBRARIES 
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${OpenMP_C_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES})
else()
    set(ADDON_LIBRARIES 
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})
endif()

set(ADDON_DEBUG_LIBRARIES
    ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG} 
    ${GSL_LIBRARY_DEBUG}
)

set(ADDON_OPTIMIZED_LIBRARIES
    ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
    ${GSL_LIBRARY}
)

add_library(hgwrbml hlmgwr.cpp)
target_link_libraries(hgwrbml ${ADDON_LIBRARIES} debug ${ADDON_DEBUG_LIBRARIES} optimized ${ADDON_OPTIMIZED_LIBRARIES})
