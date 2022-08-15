cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME dolfinx_cpx)
project(${PROJECT_NAME})

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#Set the output location.
set(BIN "${CMAKE_SOURCE_DIR}/../bin")
set(EXECUTABLE_OUTPUT_PATH ${BIN} CACHE PATH "Build directory" FORCE)
set(LIBRARY_OUTPUT_PATH ${BIN} CACHE PATH "Build directory" FORCE)

# Get DOLFIN configuration data (DOLFINXConfig.cmake must be in
# DOLFINX_CMAKE_CONFIG_PATH)
find_package(DOLFINX REQUIRED)
find_package(Basix REQUIRED)

if (PETSC_SCALAR_COMPLEX EQUAL 1)
  set(SCALAR_TYPE "--scalar_type=double _Complex")
endif()

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Compile UFL files
set(ufl_files Structural TransientThermal)
set(PYTHONPATH "$ENV{PYTHONPATH}:${CMAKE_CURRENT_SOURCE_DIR}")
foreach(ufl_file IN LISTS ufl_files)
  add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${ufl_file}.c
  COMMAND PYTHONPATH=${PYTHONPATH} python3 -m ffcx ${SCALAR_TYPE} ${CMAKE_CURRENT_SOURCE_DIR}/${ufl_file}.py
  VERBATIM DEPENDS ${ufl_file}.py problem_specification.py COMMENT "Run ffcx on ${ufl_file}.py"
)
endforeach()

# Executable
set(c_files)
foreach(ufl_file IN LISTS ufl_files)
  list(APPEND c_files ${ufl_file}.c)
endforeach()

add_library(${PROJECT_NAME} main.cpp ${c_files})

# Find Boost program_options
if(DEFINED ENV{BOOST_ROOT} OR DEFINED BOOST_ROOT)
  set(Boost_NO_SYSTEM_PATHS on)
endif()
set(Boost_USE_MULTITHREADED $ENV{BOOST_USE_MULTITHREADED})
set(Boost_VERBOSE TRUE)
find_package(Boost 1.70 REQUIRED program_options)

# Target libraries
#target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} dolfinx Boost::program_options)
target_include_directories(${PROJECT_NAME} PRIVATE ../src_op ../src ../src/Kernels)

#install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
