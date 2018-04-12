cmake_minimum_required(VERSION 3.9)
project(intersections)
set(CMAKE_CXX_STANDARD 17)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

add_executable(figures_cov catch.cpp tests_figures.cpp figures.cpp)
add_executable(figures_test catch.cpp tests_figures.cpp figures.cpp)

target_compile_options(figures_cov PRIVATE -g3 -fprofile-arcs -ftest-coverage )
set_target_properties(figures_cov PROPERTIES LINK_FLAGS "${LINK_FLAGS} -lgcov")
target_compile_options(figures_test PRIVATE -flto)