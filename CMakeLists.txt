CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mandos)

# Submodules
find_package(Git)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/eigen/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

file(GLOB SOURCE_FILES
  src/*.cpp
  src/simulation_generators/*.cpp
)
file(GLOB EDITOR_SOURCE_FILES
  src/editor/*.cpp
)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/raylib/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/tetgen/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/raylib_imgui/)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

add_library(libmandos STATIC ${SOURCE_FILES})
set_target_properties(libmandos PROPERTIES OUTPUT_NAME mandos)

add_executable(${PROJECT_NAME} src/standalone/main.cpp)
add_executable(mandos_gui ${EDITOR_SOURCE_FILES})

target_include_directories(libmandos PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/raylib_imgui/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyobj/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/raylib/src/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/tetgen/
)

target_link_libraries(libmandos PUBLIC
  raylib
  raylib_imgui
  tetgen
)

target_include_directories(${PROJECT_NAME} PUBLIC)
target_link_libraries(${PROJECT_NAME} PUBLIC libmandos)
target_include_directories(mandos_gui PUBLIC)
target_link_libraries(mandos_gui PUBLIC libmandos)
