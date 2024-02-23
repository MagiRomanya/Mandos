project(raylib_imgui)

set(public_headers
  raylib_imgui.hpp
    )

file(GLOB SOURCE_FILES
  ../imgui/*.cpp
  ../imgui/imgui_config.h
  ../imgui/backends/imgui_impl_glfw.cpp
  ../imgui/backends/imgui_impl_opengl3.cpp
  # ImGuizmo
  ../ImGuizmo/ImGuizmo.cpp
  *.cpp
)

message(${SOURCE_FILES})

add_library(raylib_imgui STATIC ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../imgui/
  ${CMAKE_CURRENT_SOURCE_DIR}/../ImGuizmo/
  ${CMAKE_CURRENT_SOURCE_DIR}/../raylib/src/
)
target_include_directories(${PROJECT_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../raylib/src/external/glfw/include
)
set_target_properties(${PROJECT_NAME} PROPERTIES
                      PUBLIC_HEADER "${public_headers}")
