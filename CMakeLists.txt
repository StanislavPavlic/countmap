cmake_minimum_required(VERSION 3.2)
project(countmap)

set(countmap_VERSION_MAJOR 1)
set(countmap_VERSION_MINOR 0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXEC_FILES src/countmap.cc src/brown_minimizers.cc src/ksw2_extz2_sse.c
               src/index.cc src/map.cc src/radixsort.cc src/region.cc
               src/revcomp.cc src/sam.cc)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/countmap.hpp.in"
  "${PROJECT_BINARY_DIR}/countmap.hpp"
)

add_library(brown_minimizers STATIC
    src/brown_minimizers.cc)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
add_subdirectory(vendor/edlib EXCLUDE_FROM_ALL)

add_executable(countmap ${EXEC_FILES})
target_link_libraries(countmap brown_minimizers bioparser thread_pool edlib pthread)
