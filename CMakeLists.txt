cmake_minimum_required(VERSION 3.9)
project(intersections)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")

set(SRC figures.cpp)
set(TEST_SRC ${SRC} catch.cpp tests_figures.cpp)

add_executable(figures_test ${TEST_SRC})
