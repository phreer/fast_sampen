cmake_minimum_required(VERSION 3.5)
project(
    kdtree_mddc
    VERSION 1.0
    LANGUAGES CXX)

message(STATUS "BUild Type: " ${CMAKE_BUILD_TYPE})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)


add_subdirectory(src)