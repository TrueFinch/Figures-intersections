cmake_minimum_required(VERSION 3.9)
project(intersections)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

add_executable(figures_cov catch.cpp figures.cpp figures.h tests.cpp)
add_executable(figures_test catch.cpp tests.cpp figures.cpp figures.h)

target_compile_options(figures_cov PRIVATE -g -O0 -Wall -fprofile-arcs -ftest-coverage)
set_target_properties(figures_cov PROPERTIES LINK_FLAGS "${LINK_FLAGS} -fprofile-arcs")