cmake_minimum_required(VERSION 3.9)
project(intersections)

#set(GCC_COMPILE_FLAGS "-Wall -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")
set(CMAKE_CXX_STANDARD 14)

#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)

add_executable(figures_cov catch.cpp tests.cpp figures.cpp figures.cpp figures.h tests.cpp)
add_executable(figures_test catch.cpp tests.cpp figures.cpp figures.h)

target_compile_options(figures_cov PRIVATE -g3 -fprofile-arcs -ftest-coverage)
set_target_properties(figures_cov PROPERTIES LINK_FLAGS "${LINK_FLAGS} -lgcov")
#target_compile_options(figures_test PRIVATE -flto)