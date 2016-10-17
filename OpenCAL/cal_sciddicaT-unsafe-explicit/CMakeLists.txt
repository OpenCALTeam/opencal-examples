cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_package(OpenCAL 1.1 COMPONENTS  OpenCAL  REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${OPENCAL_INCLUDE_DIR})

file(GLOB SRC
    "include/*.h"
    "source/*.c"
)

add_executable(cal_sciddicaT-unsafe-explicit ${SRC})
target_link_libraries(cal_sciddicaT-unsafe-explicit ${OPENCAL_LIBRARIES})
