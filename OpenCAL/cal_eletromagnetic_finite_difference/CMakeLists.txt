# cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#
# find_package(OpenCAL 1.1 COMPONENTS  OpenCAL OpenCAL-GL REQUIRED)
# find_package(Threads REQUIRED)
# find_package(OpenGL  REQUIRED)
# find_package(GLUT REQUIRED)
#
# if(NOT MSVC)
#   find_library(M_LIB m)
# else()
#   # not needed on windows
#   set(M_LIB "")
# endif()
#
# # FIND_PACKAGE(OpenMP)
# # IF(OPENMP_FOUND)
# #     SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# #     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# # ENDIF()
#
# include_directories(
# 	${CMAKE_CURRENT_SOURCE_DIR}/include
# 	${OPENGL_INCLUDE_DIRS}
# 	${GLUT_INCLUDE_DIR}
# 	${OPENCAL_INCLUDE_DIR}
# )
#
# link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
#
# file(GLOB SRC
#     "include/*.h"
#     "source/*.cpp"
#     "source/*.c"
# )
#
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
# add_executable(cal_eletromagnetic_finite_difference ${SRC})

# target_link_libraries(cal_eletromagnetic_finite_difference  ${OPENCAL_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
