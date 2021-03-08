cmake_minimum_required (VERSION 3.14)

project (lgrid CXX)

add_executable (lgrid)
add_dependencies (lgrid stell)

target_link_libraries (lgrid stell)

add_subdirectory (Sources)
