cmake_minimum_required(VERSION 3.22.0)
project(dmpermlib VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

# Ensure Release build type is set
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add dependecies
add_subdirectory(src/extern/CSparse EXCLUDE_FROM_ALL)

# nanobind
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

# Find the NumPy include directory
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import numpy as np; print(np.get_include())"
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NUMPY_DIR)

# Add the extension module to the build
nanobind_add_module(dmpermlib src/src_extension/dmpermlib.cc src/src_extension/SparseMatrix.cc src/src_extension/StructuralAnalysisModel.cc src/src_extension/MSOAlg.cc)
target_link_libraries(dmpermlib PRIVATE csparse)
target_include_directories(dmpermlib PRIVATE ${NUMPY_DIR})
target_compile_features(dmpermlib PRIVATE cxx_std_20)
if (NOT WIN32)
  target_compile_options(dmpermlib PRIVATE -Wall -Wextra -Wno-nested-anon-types -Wno-zero-length-array)
endif()

# Copy the extension module to the src directory
# add_custom_command(TARGET dmpermlib POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:dmpermlib> ${PROJECT_SOURCE_DIR}/src/faultdiagnosistoolbox)

# install the generated extension module in the toolbox directory
install(TARGETS dmpermlib LIBRARY DESTINATION faultdiagnosistoolbox)
