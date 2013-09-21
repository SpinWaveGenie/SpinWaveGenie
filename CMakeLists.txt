project(spin_wave_genie)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE")

find_package(Eigen3 REQUIRED)
find_package(Boost 1.5 REQUIRED COMPONENTS system thread program_options)

add_subdirectory(src)
