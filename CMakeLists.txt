project(spin_wave_genie)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

cmake_minimum_required(VERSION 2.8)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

ENABLE_TESTING()


macro(remove_cxx_flag flag)
  string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endmacro()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#find_package(Log4cplus)
#if(NOT ${LOG4CPLUS_FOUND})
#  include(ExternalProject)
#  ExternalProject_Add(log4cplus
#    GIT_REPOSITORY https://github.com/log4cplus/log4cplus.git
#    GIT_TAG REL_1_1_2
#    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/log4cplus
#    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/log4cplus "-DCMAKE_CXX_FLAGS:STRING=-std=c++11 -stdlib=libc++" "-DCMAKE_EXE_LINKER_FLAGS:STRING=-std=c++11 -stdlib=libc++" -DBUILD_SHARED_LIBS:BOOL=FALSE 
#)
#  set(LOG4CPLUS_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/log4cplus/include)
#  set(LOG4CPLUS_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/log4cplus/lib/liblog4cplus.dylib)
#endif(NOT ${LOG4CPLUS_FOUND})

#message("sfml_PREFIX='${LOG4CPLUS_LIBRARIES}'")

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
find_package(TBB REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.5 REQUIRED COMPONENTS system thread program_options unit_test_framework)

add_subdirectory(src)
add_subdirectory(test)
 

