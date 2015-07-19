cmake_minimum_required(VERSION 3.2)
project(regularization)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(gmock)


set(SOURCE_FILES main.cpp)

add_executable(regularization ${SOURCE_FILES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(PopulationCore)

target_link_libraries(regularization populationCore)
