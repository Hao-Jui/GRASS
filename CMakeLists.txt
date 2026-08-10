cmake_minimum_required(VERSION 3.21)
project(GRASS Fortran)

# Compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG
  "-O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -fcheck=all -Wconversion -finit-real=nan")

# Sanitizer build type (ASan + UBSan)
set(CMAKE_Fortran_FLAGS_SANITIZER
  "-O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -fcheck=all -finit-real=nan -fsanitize=address,undefined")
set(CMAKE_EXE_LINKER_FLAGS_SANITIZER "-fsanitize=address,undefined")

# External dependencies
if(APPLE)
  set(BLA_VENDOR Apple)
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

# Module output directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Source files (all except main.f90)
file(GLOB CORE_SOURCES src/core/*.f90)
file(GLOB TOOL_SOURCES src/tool/*.f90)
file(GLOB THEORY_SOURCES src/theory/*.f90)
file(GLOB BH_TOROID_SOURCES src/bh_toroid/*.f90)
file(GLOB UTILITIES_SOURCES src/utilities/*.f90)
set(LIB_SOURCES ${CORE_SOURCES} ${TOOL_SOURCES} ${THEORY_SOURCES} ${BH_TOROID_SOURCES}
  ${UTILITIES_SOURCES}
  src/para_panel.f90)

# Library target
add_library(grass_lib STATIC ${LIB_SOURCES})
set_target_properties(grass_lib PROPERTIES
  OUTPUT_NAME grass
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
target_include_directories(grass_lib PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${FFTW3_INCLUDE_DIRS})
target_link_directories(grass_lib PUBLIC ${FFTW3_LIBRARY_DIRS})
target_link_libraries(grass_lib PUBLIC
  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${FFTW3_LIBRARIES})

# Main executable
add_executable(grass src/main.f90)
target_link_libraries(grass grass_lib)

# Platform-conditional stack size
if(APPLE)
  target_link_options(grass PRIVATE "-Wl,-stack_size,0x4000000")
endif()

# Install targets
include(GNUInstallDirs)
install(TARGETS grass RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS grass_lib ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Tests
enable_testing()
add_subdirectory(tests)
