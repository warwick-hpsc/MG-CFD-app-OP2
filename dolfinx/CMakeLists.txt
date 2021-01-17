Project(CPX)
cmake_minimum_required(VERSION 3.10)

#Set location of our FindFOO.cmake modules
#find_path(CPX_CMAKE_DIR modules PATHS "${CMAKE_SOURCE_DIR}/cmake" NO_DEFAULT_PATH)
#set(CMAKE_MODULE_PATH ${CPX_CMAKE_DIR}/modules CACHE FILEPATH
#     "search path for scripts and modules")

#Set the output location.
set(BIN "${CMAKE_SOURCE_DIR}/../bin")
set(EXECUTABLE_OUTPUT_PATH ${BIN} CACHE PATH "Build directory" FORCE)
set(LIBRARY_OUTPUT_PATH ${BIN} CACHE PATH "Build directory" FORCE)

# Get DOLFINX configuration data (DOLFINXConfig.cmake must be in
# DOLFINX_CMAKE_CONFIG_PATH)
find_package(DOLFINX REQUIRED)
find_package(HDF5 REQUIRED)

# Executable
add_library(dolfinx_cpx STATIC "main.cpp" "poisson.c")
#add_executable(${PROJECT_NAME} entry.cpp poisson.c main.cpp)

# Set C++17 standard
target_compile_features(dolfinx_cpx PRIVATE cxx_std_17)

if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic") #-Werror")
    set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS})
endif()

# Target libraries
target_link_libraries(dolfinx_cpx dolfinx)
target_include_directories(dolfinx_cpx PRIVATE "../src" "../src/Kernels")
# Do not throw error for 'multi-line comments' (these are typical in
# rst which includes LaTeX)
#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-Wno-comment" HAVE_NO_MULTLINE)
#set_source_files_properties(entry.cpp PROPERTIES COMPILE_FLAGS "$<$<BOOL:${HAVE_NO_MULTLINE}>:-Wno-comment -Wall -Wextra -pedantic -Werror>")
