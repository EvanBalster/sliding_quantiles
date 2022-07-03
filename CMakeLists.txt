cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
# Currently need C++17 for std::size, otherwise C++14 works
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_language(CXX)
enable_language(C)

find_package (Threads)

# set the project name
project(SlidingQuantiles)

file(GLOB QUERN_HEADERS include/quern/*.hpp)

# add the executable
add_executable(SlidingQuantiles test/main.cpp ${QUERN_HEADERS})

target_include_directories(SlidingQuantiles PUBLIC "include")
