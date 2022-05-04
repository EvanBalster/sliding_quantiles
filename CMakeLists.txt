cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(SlidingQuantiles)

# add the executable
add_executable(SlidingQuantiles test/main.cpp include/histogram_quantiles.hpp)

target_include_directories(SlidingQuantiles PUBLIC "include")
