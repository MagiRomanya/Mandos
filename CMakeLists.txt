cmake_minimum_required(VERSION 3.9)
project(mandos LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# CXX Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

# Make visibility hidden by default, as in Windows.
# Its the recommended way to separate private and public API 
# https://www.youtube.com/watch?v=m0DwB4OvDXk
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

include(GenerateExportHeader)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
# find_package(pybind11 REQUIRED)
find_package(spdlog REQUIRED)
find_package(Tracy REQUIRED)

# Support interface target for adding project wise options
add_custom_target(project_options)

# Link time optimization support
include(CheckIPOSupported)
check_ipo_supported(RESULT link_time_optimization_support)

# Enable LTO for the project
if(link_time_optimization_support)
  set_target_properties(project_options PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Set output directories for libraries, archives and binaries
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(src)
# add_subdirectory(bindings)

if (BUILD_TESTING)
  find_package(Catch2 REQUIRED)

endif()