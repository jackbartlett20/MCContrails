cmake_minimum_required(VERSION 3.10)
project(MCContrails LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Release (default) or Debug" FORCE)
endif()

# Check for libraries
find_package(yaml-cpp REQUIRED)
find_package(OpenMP REQUIRED)

file(GLOB_RECURSE SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(MCContrails ${SOURCE_FILES})

# Add include directories
target_include_directories(MCContrails PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Include libraries
target_link_libraries(MCContrails PRIVATE yaml-cpp::yaml-cpp)
target_link_libraries(MCContrails PUBLIC OpenMP::OpenMP_CXX)

# Compile features and options
target_compile_features(MCContrails PRIVATE cxx_std_20)
target_compile_options(MCContrails PRIVATE -Wall -Wextra)
target_compile_options(MCContrails PUBLIC -fopenmp)
# Flags for Debug build configuration
target_compile_options(MCContrails PRIVATE $<$<CONFIG:Debug>:-O0> $<$<CONFIG:Debug>:-g>)
# Flags for Release build configuration (default)
target_compile_options(MCContrails PRIVATE $<$<CONFIG:Release>:-O3>)