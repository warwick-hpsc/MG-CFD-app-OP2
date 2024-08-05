cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME refine_mesh)
project(${PROJECT_NAME})

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(Boost_USE_STATIC_LIBS ON)
# Get DOLFIN configuration data (DOLFINXConfig.cmake must be in
# DOLFINX_CMAKE_CONFIG_PATH)
find_package(DOLFINX REQUIRED)
find_package(Basix REQUIRED)

if (PETSC_SCALAR_COMPLEX EQUAL 1)
  set(SCALAR_TYPE "--scalar_type=double _Complex")
endif()

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(${PROJECT_NAME} main.cpp)

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
#target_include_directories(${PROJECT_NAME} PRIVATE ../src_op ../src ../src/Kernels)
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
#install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
