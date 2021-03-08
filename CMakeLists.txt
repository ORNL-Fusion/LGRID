cmake_minimum_required (VERSION 3.14)

project (lgrid CXX)

set (CMAKE_CXX_STANDARD 11)

add_executable (lgrid)
add_dependencies (lgrid stell)

target_link_libraries (lgrid stell)

add_subdirectory (Sources)
