cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

find_package(OpenCAL 1.1 COMPONENTS  OpenCAL OpenCAL-CL  REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

find_package( OpenCL REQUIRED )

if(NOT MSVC)
  find_library(M_LIB m)
else()
  # not needed on windows
  set(M_LIB "")
endif()
include_directories( ${OpenCL_INCLUDE_DIRS} )
link_directories(${OpenCL_LIBRARIES})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${OPENCAL_INCLUDE_DIR})

file(GLOB SRC
    "include/*.h"
    "source/*.c"
)

add_executable(calcl_life ${SRC})
target_link_libraries(calcl_life ${OPENCAL_LIBRARIES} ${M_LIB} ${OpenCL_LIBRARY})
