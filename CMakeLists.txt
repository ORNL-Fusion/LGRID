cmake_minimum_required (VERSION 3.14)

project (lgrid CXX)

add_executable (lgrid)

target_link_libraries (lgrid PUBLIC stell)
target_compile_features (lgrid PUBLIC cxx_std_11)

add_subdirectory (Sources)
