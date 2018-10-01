cmake_minimum_required(VERSION 2.8)

project(AHP-Matrix-Project)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

file(GLOB_RECURSE common ${CMAKE_SOURCE_DIR}/sources/common/*.*)
file(GLOB_RECURSE processor ${CMAKE_SOURCE_DIR}/sources/projects/higherDimRandom/*.*
                            ${CMAKE_SOURCE_DIR}/sources/projects/higherDimRandom/routines/*.*
                            ${CMAKE_SOURCE_DIR}/sources/projects/higherDimRandom/counters/*.*)

file(GLOB_RECURSE processorMultithreaded ${CMAKE_SOURCE_DIR}/sources/projects/higherDimRandomMultithreaded/*.*
                                         ${CMAKE_SOURCE_DIR}/sources/projects/higherDimRandomMultithreaded/routines/*.*
                                         ${CMAKE_SOURCE_DIR}/sources/projects/higherDimRandomMultithreaded/counters/*.*)

#add_subdirectory(include/lemon)

include_directories(include
                    sources/common
                    sources/projects
                    sources/projects/higherDimRandom/counters
                    sources/projects/higherDimRandomMultithreaded/counters
                    sources/projects/higherDimRandomMultithreaded/routines)

#Flags and standard setting globally
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -o2 -Wall -Wextra -Werror") #-lglpk -lemon

add_executable(main ${common} sources/projects/main/main.cpp)
add_executable(dimensionFive ${common} sources/projects/dimensionFive/Gen5x5Matrices.cpp)
add_executable(dimensionFiveThreaded ${common} sources/projects/dimensionFive/Multithreaded.cpp)
add_executable(dimensionFiveThreaded2 ${common} sources/projects/dimensionFive/Multithreaded2.cpp)
add_executable(dimensionFiveThreaded3 ${common} sources/projects/dimensionFive/Multithreaded3.cpp)
add_executable(higherDimRandom ${common} ${processor} sources/projects/higherDimRandom/higherDimRandom.cpp)
add_executable(higherDimRandomMultithreaded ${common} ${processorMultithreaded} sources/projects/higherDimRandomMultithreaded/higherDimRandomMultithreaded.cpp)
add_executable(histogram ${common} sources/projects/histogram/Histogram.cpp)
add_executable(lpTester ${common} sources/projects/test/LPTester.cpp)

TARGET_LINK_LIBRARIES(main emon glpk)
TARGET_LINK_LIBRARIES(dimensionFive emon glpk)
TARGET_LINK_LIBRARIES(dimensionFiveThreaded emon glpk pthread)
TARGET_LINK_LIBRARIES(dimensionFiveThreaded2 emon glpk pthread)
TARGET_LINK_LIBRARIES(dimensionFiveThreaded3 emon glpk pthread)
TARGET_LINK_LIBRARIES(higherDimRandom emon glpk)
TARGET_LINK_LIBRARIES(higherDimRandomMultithreaded emon glpk pthread)
TARGET_LINK_LIBRARIES(histogram emon glpk)
TARGET_LINK_LIBRARIES(lpTester emon glpk)
