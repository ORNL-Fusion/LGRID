target_sources(lgrid
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/main.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vector3d.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vector3d.hpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vertex.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vertex.hpp>
)
