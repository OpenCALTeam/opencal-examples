cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

find_package(OpenCAL 1.0 COMPONENTS  OpenCAL OpenCAL-CL  REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

find_package( OpenCL REQUIRED )
include_directories( ${OpenCL_INCLUDE_DIRS}
                     ${OPENCAL_INCLUDE_DIR}
                   )
link_directories(${OpenCL_LIBRARIES})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB SRC
    "include/*.h"
    "source/*.cpp"
)

add_executable(calcl_mbusu ${SRC})
target_link_libraries(calcl_mbusu ${OPENCAL_LIBRARIES} ${OpenCL_LIBRARY})
