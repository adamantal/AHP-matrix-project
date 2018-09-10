cmake_minimum_required(VERSION 2.8)

project(AHP-Matrix-Project)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

file(GLOB_RECURSE common ${CMAKE_SOURCE_DIR}/sources/common/*.cpp)

include_directories(include)
include_directories(sources/common)
include_directories(sources/projects)

#Flags and standard setting globally
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -o3 -Wall -Wextra -Werror") #-lglpk -lemon

add_executable(main ${common} sources/projects/main/main.cpp)
add_executable(dimensionFive ${common} sources/projects/dimensionFive/Gen5x5Matrices.cpp)
add_executable(histogram ${common} sources/projects/histogram/Histogram.cpp)
add_executable(lpTester ${common} sources/projects/test/LPTester.cpp)

TARGET_LINK_LIBRARIES(main emon glpk)
TARGET_LINK_LIBRARIES(dimensionFive emon glpk)
TARGET_LINK_LIBRARIES(histogram emon glpk)
TARGET_LINK_LIBRARIES(lpTester emon glpk)
