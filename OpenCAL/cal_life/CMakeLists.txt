# opencal/OpenCAL/examples/cal_life/CMakeLists.txt
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_package(OpenCAL 1.0 COMPONENTS OpenCAL REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${OPENCAL_INCLUDE_DIR}
                    )

file(GLOB SRC
    "include/*.h"
    "source/*.c"
)

add_executable(cal_life ${SRC})
target_link_libraries(cal_life ${OPENCAL_LIBRARIES})
