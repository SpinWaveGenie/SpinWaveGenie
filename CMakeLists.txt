project(spin_wave_genie)

cmake_minimum_required(VERSION 2.8)

macro(remove_cxx_flag flag)
  string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endmacro()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE -msse4")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -DEIGEN_DONT_PARALLELIZE -msse4 -O4")
remove_cxx_flag("-O3")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++11 -stdlib=libc++")
find_package(Eigen3 REQUIRED)
find_package(Boost 1.5 REQUIRED COMPONENTS system thread program_options)

add_subdirectory(src)
