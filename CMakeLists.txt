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

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

macro(remove_cxx_flag flag)
  string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endmacro()

if(CMAKE_COMPILER_IS_GNUCXX)
    message("compiler detected as: GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE -std=c++0x -march=native")
    #remove_cxx_flag("-O3")
elseif(APPLE)
    message("compiler detected as: LLVM compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -DEIGEN_DONT_PARALLELIZE -msse4 -ffast-math")
endif(CMAKE_COMPILER_IS_GNUCXX)

#find_package(Log4cplus)
find_package(NLOPT REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(Nexus REQUIRED)
#set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.5 REQUIRED COMPONENTS system thread program_options unit_test_framework)

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

add_subdirectory(src)
add_subdirectory(test)
 

