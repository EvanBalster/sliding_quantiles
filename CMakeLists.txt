cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(SlidingQuantiles)

file(GLOB QUERN_HEADERS include/quern/*.hpp)

# add the executable
add_executable(SlidingQuantiles test/main.cpp ${QUERN_HEADERS})

target_include_directories(SlidingQuantiles PUBLIC "include")
