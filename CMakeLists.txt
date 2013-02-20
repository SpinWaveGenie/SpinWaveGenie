project(spin_wave_genie)

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)

add_subdirectory(src)
