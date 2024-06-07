CMAKE_MINIMUM_REQUIRED(VERSION 3.9)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mandos)

option(USE_VISUALIZER "Compile the built in simulation visualizer alongside the simulation" ON)
option(USE_PYTHON_LIBRARY "Generate the python bindings for using the simulation library" ON)

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

# Configure file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/mandos_config.hpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mandos_config.hpp
)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/eigen/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # Necessary for pybind11 dynamic lib
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/tetgen/)

set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_BUILD_TYPE Release)
# add_compile_options(-fno-omit-frame-pointer) # for profiling
find_package(OpenMP)

file(GLOB SOURCE_FILES
  src/*.cpp
  src/simulation_generators/*.cpp
  src/optimization/*.cpp
)

if (USE_VISUALIZER)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/raylib/)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/raylib_imgui/)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/nativefiledialog-extended/)
  file(GLOB RENDERING_SOURCE_FILES
    src/rendering/*.cpp
  )
endif()

add_library(libmandos STATIC ${SOURCE_FILES} ${RENDERING_SOURCE_FILES})
set_target_properties(libmandos PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})

# Link time optimization support
include(CheckIPOSupported)
check_ipo_supported(RESULT link_time_optimization_support)

if(link_time_optimization_support)
  set_target_properties(libmandos PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

target_include_directories(libmandos PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyobj/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/tetgen/
  ${CMAKE_CURRENT_SOURCE_DIR}/external/TriangleMeshDistance/TriangleMeshDistance/include/
)

target_link_libraries(libmandos
  tetgen
  # profiler
  # OpenMP::OpenMP_CXX
)


if (USE_VISUALIZER)
  target_include_directories(libmandos PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/
    ${CMAKE_CURRENT_SOURCE_DIR}/external/raylib_imgui/
    ${CMAKE_CURRENT_SOURCE_DIR}/external/raylib/src/
    ${CMAKE_CURRENT_SOURCE_DIR}/external/nativefiledialog-extended/src/include/
  )

  target_link_libraries(libmandos
    raylib
    raylib_imgui
    nfd
  )

# Build the main executables if the visualizer is shipped with the library
  add_executable(${PROJECT_NAME} src/standalone/main.cpp)
  add_executable(mandos2 src/standalone/main2.cpp)

  target_link_libraries(${PROJECT_NAME} PUBLIC libmandos)
  target_link_libraries(mandos2 PUBLIC libmandos)
endif()

if (USE_PYTHON_LIBRARY)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11/)

  # Creating the pyhon library
  find_package(Python 3.6 COMPONENTS Interpreter Development REQUIRED)
  pybind11_add_module(pymandos ${CMAKE_CURRENT_SOURCE_DIR}/src/python/pybinding.cpp)
  target_link_libraries(pymandos PUBLIC libmandos ${PYTHON_LIBRARIES})
  include_directories(pymandos PUBLIC ${PYTHON_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11/include/
  )
endif()
