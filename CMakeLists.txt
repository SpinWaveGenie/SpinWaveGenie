project(spin_wave_genie)

cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(SPINWAVEGENIE_VERSION_MAJOR "0")
set(SPINWAVEGENIE_VERSION_MINOR "1")
set(SPINWAVEGENIE_VERSION_PATCH "0")
set(SPINWAVEGENIE_VERSION "${SPINWAVEGENIE_VERSION_MAJOR}.${SPINWAVEGENIE_VERSION_MINOR}.${SPINWAVEGENIE_VERSION_PATCH}")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

if(CMAKE_COMPILER_IS_GNUCXX)
    message("compiler detected as: GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE -std=c++0x -march=native" CACHE STRING "compile flags" FORCE)
elseif(APPLE)
    message("Setting up Xcode for C++11 with libc++.")
    if(CMAKE_GENERATOR STREQUAL Xcode)
        set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++0x")
        set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -DEIGEN_DONT_PARALLELIZE -msse4" CACHE STRING "compile flags" FORCE)
endif(CMAKE_COMPILER_IS_GNUCXX)

find_package(NLOPT REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(Nexus REQUIRED)

string(FIND ${CMAKE_SYSTEM} "el6" isRHEL6)
if (isRHEL6)
    set(Boost_NO_BOOST_CMAKE ON)
endif(isRHEL6)
find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)

option(BUILD_EXAMPLES "Build SpinWaveGenieExamples" OFF)
option(BUILD_TESTING "Build unit tests" OFF)
option(ENABLE_DOXYGEN "Enables building of documentation with doxygen" OFF)

if(BUILD_TESTING)
    enable_testing()
endif(BUILD_TESTING)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(ENABLE_DOXYGEN AND DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(ENABLE_DOXYGEN AND DOXYGEN_FOUND)

add_subdirectory(libSpinWaveGenie)

