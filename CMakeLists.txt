cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

option(BUILD_OPENCAL        "Build the OpenCAL serial version"                                      ON)
option(BUILD_OPENCAL_OMP    "Build the OpenCAL-OMP OpenMP parallel version (OpenMP required)"       ON)
option(BUILD_OPENCAL_CL     "Build the OpenCAL-CL OpenCL parallel version (OpenCL required)"        ON)
option(BUILD_OPENCAL_GL     "Build the OpenCAL-GL visualization library (OpenGL and GLUT required)" ON)
option(BUILD_OPENCAL_++     "Build the OpenCAL-GL visualization library (OpenGL and GLUT required)" ON)


if(BUILD_OPENCAL)
  add_subdirectory(OpenCAL)
endif(BUILD_OPENCAL)

if(BUILD_OPENCAL_OMP)
add_subdirectory(OpenCAL-OMP)
endif(BUILD_OPENCAL_OMP)

if(BUILD_OPENCAL_GL)
add_subdirectory(OpenCAL-GL)
endif(BUILD_OPENCAL_GL)

if(BUILD_OPENCAL_CL)
add_subdirectory(OpenCAL-CL)
endif(BUILD_OPENCAL_CL)

if(BUILD_OPENCAL_++)
add_subdirectory(OpenCAL++)
endif(BUILD_OPENCAL_++)
