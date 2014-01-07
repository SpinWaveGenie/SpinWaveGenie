project(spin_wave_genie)

cmake_minimum_required(VERSION 2.8)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

ENABLE_TESTING()

macro(remove_cxx_flag flag)
  string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endmacro()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if(CMAKE_COMPILER_IS_GNUCXX)
    message("compiler detected as: GNU compiler")
    set(EIGEN3_INCLUDE_DIR "/SNS/users/svh/software/eigen3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE -std=c++0x -march=native")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++0x -march=native") 
    set(EIGEN3_INCLUDE_DIR "/SNS/users/svh/software/eigen3")
elseif(APPLE)
    message("compiler detected as: LLVM compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -DEIGEN_DONT_PARALLELIZE -msse4")
    #remove_cxx_flag("-O3")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++11 -stdlib=libc++ -msse4")
endif(CMAKE_COMPILER_IS_GNUCXX)
find_package(Eigen3 REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.5 REQUIRED COMPONENTS system thread program_options unit_test_framework)

add_subdirectory(src)
add_subdirectory(test)
 
