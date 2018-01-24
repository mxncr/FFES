PROJECT(geogram_basic C CXX)
cmake_minimum_required(VERSION 3.0)

message(STATUS ${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/..)
file(GLOB LIB_SOURCES "basic/*.cpp" "api/*.cpp")
file(GLOB LIB_HEADERS "basic/*.h" "api/.h")
add_library(geogram_basic STATIC ${LIB_SOURCES} ${LIB_HEADERS})
