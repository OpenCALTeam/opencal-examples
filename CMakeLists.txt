cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(OPENCAL_MAJOR_VERSION 1)
set(OPENCAL_MINOR_VERSION 0)
set(OPENCAL_VERSION
  ${OPENCAL_MAJOR_VERSION}.${OPENCAL_MINOR_VERSION})

project(OpenCAL-examples-${OPENCAL_VERSION}
		LANGUAGES C CXX)

option(BUILD_OPENCAL_SERIAL "Build the OpenCAL serial version"                                      ON )
option(BUILD_OPENCAL_OMP    "Build the OpenCAL-OMP OpenMP parallel version (OpenMP required)"       ON )
option(BUILD_OPENCAL_OMP_PARALLEL "Controls if OpenCAL-OMP compiled agaist libomp. If off just one processor is used."  ON )
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
