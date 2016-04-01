# opencal/OpenCAL-CL/examples/calcl_gl_mod2CA3D/CMakeLists.txt
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_package(OpenCAL 1.0 COMPONENTS  OpenCAL OpenCAL-CL OpenCAL-GL REQUIRED)
find_package(OpenCL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Threads REQUIRED)
find_library(M_LIB m)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
  ${OPENCAL_INCLUDE_DIR}
	${OPENGL_INCLUDE_DIRS}
	${GLUT_INCLUDE_DIR}
)
include_directories( ${OpenCL_INCLUDE_DIRS} )

link_directories(${OpenCL_LIBRARIES})

file(GLOB SRC
    "include/*.h"
    "source/*.cpp"
    "source/*.c"
)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(calcl_gl_mod2CA3D ${SRC})
target_link_libraries(calcl_gl_mod2CA3D ${OPENCAL_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} -lpthread ${OpenCL_LIBRARY})
