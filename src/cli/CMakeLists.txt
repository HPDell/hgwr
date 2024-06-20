cmake_minimum_required(VERSION 3.12.0)

include_directories(../include)

find_package(Boost COMPONENTS program_options)
if(Boost_FOUND)
    message(STATUS "Boost found")
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost Library Cannot Found!")
endif()

if(DEFINED Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG)
    set(ADDON_DEBUG_LIBRARIES
        ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG} 
    )
else()
    set(ADDON_DEBUG_LIBRARIES
        ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE} 
    )
endif()
set(ADDON_OPTIMIZED_LIBRARIES
    ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
)


if(NOT MSVC)
    add_definitions(-DDBL_MAX=__DBL_MAX__)
endif()

add_executable(hgwr main.cpp)
target_link_libraries(hgwr PRIVATE hgwrbml debug ${ADDON_DEBUG_LIBRARIES} optimized ${ADDON_OPTIMIZED_LIBRARIES})