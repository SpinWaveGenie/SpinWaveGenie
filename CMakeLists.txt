project(spin_wave_genie)

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

find_package(Eigen3 REQUIRED)
find_package(Boost 1.5 REQUIRED COMPONENTS system thread program_options)

add_subdirectory(src)
