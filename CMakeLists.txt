# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(rcbkdipole CXX C)

# Optimization on (for debug flags, change to Debug)
set(CMAKE_BUILD_TYPE Release)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

add_subdirectory(src)

set_target_properties(amplitude PROPERTIES
	 LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
	ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

set_target_properties(dipole PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
