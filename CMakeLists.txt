cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

option(BUILD_OPENCAL_SERIAL "Build the OpenCAL serial version"                                      ON )
option(BUILD_OPENCAL_OMP    "Build the OpenCAL-OMP OpenMP parallel version (OpenMP required)"       ON )
option(BUILD_OPENCAL_CL     "Build the OpenCAL-CL OpenCL parallel version (OpenCL required)"        OFF)
option(BUILD_OPENCAL_GL     "Build the OpenCAL-GL visualization library (OpenGL and GLUT required)" OFF)
option(BUILD_OPENCAL_++     "Build the OpenCAL-GL visualization library (OpenGL and GLUT required)" OFF)


if(BUILD_OPENCAL_SERIAL)
  add_subdirectory(OpenCAL)
endif(BUILD_OPENCAL_SERIAL)

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
