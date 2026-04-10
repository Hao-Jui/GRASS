cmake_minimum_required(VERSION 3.12)
project(GRASS Fortran)

# Compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG
  "-O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -fcheck=all -Wconversion -finit-real=nan")

# External dependencies
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

# Module output directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Source files (all except main.f90)
file(GLOB CORE_SOURCES src/core/*.f90)
file(GLOB TOOL_SOURCES src/tool/*.f90)
file(GLOB THEORY_SOURCES src/theory/*.f90)
set(LIB_SOURCES ${CORE_SOURCES} ${TOOL_SOURCES} ${THEORY_SOURCES}
  src/para_panel.f90)

# Library target
add_library(grass_lib STATIC ${LIB_SOURCES})
target_include_directories(grass_lib PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${FFTW3_INCLUDE_DIRS})
target_link_directories(grass_lib PUBLIC ${FFTW3_LIBRARY_DIRS})
target_link_libraries(grass_lib PUBLIC
  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${FFTW3_LIBRARIES})

# Main executable
add_executable(grass src/main.f90)
target_link_libraries(grass grass_lib)

# Tests
enable_testing()
add_subdirectory(tests)
