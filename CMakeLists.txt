cmake_minimum_required(VERSION 2.8)

project(AHP-Matrix-Project)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

file(GLOB_RECURSE sources ${CMAKE_SOURCE_DIR}/sources/Matrices/*.cpp) #${CMAKE_SOURCE_DIR}/sources/Matrices/*.hpp

#find_package(LEMON)

include_directories(include)
include_directories(sources/Matrices)
include_directories(sources)

#Flags and standard setting globally
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -o3 -Wall -Wextra -Werror") #-lglpk -lemon

add_executable(Histogram ${sources} sources/_Histogram.cpp)

TARGET_LINK_LIBRARIES(Histogram emon)
TARGET_LINK_LIBRARIES(Histogram glpk) 
