cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_package(OpenCAL 1.1 COMPONENTS  OpenCAL OpenCAL-CL  REQUIRED)

if(NOT MSVC)
  find_library(M_LIB m)
else()
  # not needed on windows
  set(M_LIB "")
endif()

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
    "source/*.c"
)

add_executable(calcl_mod2CA3D ${SRC})
target_link_libraries(calcl_mod2CA3D ${OPENCAL_LIBRARIES} ${M_LIB} ${OpenCL_LIBRARY})
