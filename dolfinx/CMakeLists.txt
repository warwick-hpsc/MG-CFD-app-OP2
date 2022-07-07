cmake_minimum_required(VERSION 3.10)

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

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in
# DOLFIN_CMAKE_CONFIG_PATH)
#find_package(Eigen3 REQUIRED)
find_package(DOLFINX REQUIRED)
find_package(Basix REQUIRED)

set(CMAKE_BUILD_TYPE "Release")

# Compile UFL files
add_custom_command(
  OUTPUT Structural.c
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/Structural.ufl
  DEPENDS Structural.ufl
)
add_custom_command(
  OUTPUT TransientThermal.c
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/TransientThermal.ufl
  DEPENDS TransientThermal.ufl
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Executable
add_library(${PROJECT_NAME} main.cpp 
  ${CMAKE_CURRENT_BINARY_DIR}/Structural.c
  ${CMAKE_CURRENT_BINARY_DIR}/TransientThermal.c)

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
