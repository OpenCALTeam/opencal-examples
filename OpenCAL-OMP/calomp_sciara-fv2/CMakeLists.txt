cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_package(OpenCAL 1.0 COMPONENTS  OpenCAL-OMP REQUIRED)
find_library(M_LIB m)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${OPENCAL_INCLUDE_DIR})

file(GLOB SRC
    "include/*.h"
    "source/*.cpp"
)

add_executable(calomp_sciara_fv2 ${SRC})
target_link_libraries(calomp_sciara_fv2 ${M_LIB} ${OPENCAL_LIBRARIES} )
