project(SpinWaveGenie)

cmake_minimum_required (VERSION 3.1)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MACOSX_RPATH TRUE)

set(SPINWAVEGENIE_VERSION_MAJOR "0")
set(SPINWAVEGENIE_VERSION_MINOR "2")
set(SPINWAVEGENIE_VERSION_PATCH "1")
set(SPINWAVEGENIE_VERSION "${SPINWAVEGENIE_VERSION_MAJOR}.${SPINWAVEGENIE_VERSION_MINOR}.${SPINWAVEGENIE_VERSION_PATCH}")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

if(BUILD_TESTING)
    enable_testing()
endif(BUILD_TESTING)

# We probably don't want this to run on every build.
option(COVERALLS "Generate coveralls data" OFF)

if (COVERALLS)
    option(COVERALLS_TRAVISCI "Generate coveralls data for travis-ci" ON)
    option(COVERALLS_UPLOAD "Automatically upload coveralls data" ON)
    include(Coveralls)
    coveralls_turn_on_coverage()
endif()

set(CMAKE_CXX_STANDARD 14)
set(CXX_STANDARD_REQUIRED 11)
# XCode isn't picking up the standard set above.
if(CMAKE_GENERATOR STREQUAL Xcode)
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE -Wall -Wextra")
endif()

if(CMAKE_GENERATOR STREQUAL Xcode)
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
endif()
option(WITH_SANITIZERS "Enable sanitizers" OFF)
if(WITH_SANITIZERS)
  message(STATUS "enabling address & undefined behavior sanitizers")
  file(WRITE ${PROJECT_BINARY_DIR}/ubsan.txt "fun:*_Ios_Fmtflags*")
  set(SAN_FLAGS "-fno-omit-frame-pointer -fno-common -fsanitize=address,undefined -fno-sanitize-recover -fsanitize-blacklist=${PROJECT_BINARY_DIR}/ubsan.txt")
  add_compile_options(-fno-omit-frame-pointer -fno-common -fsanitize=address,undefined -fno-sanitize-recover -fsanitize-blacklist=${PROJECT_BINARY_DIR}/ubsan.txt)
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${SAN_FLAGS}" )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SAN_FLAGS}" )
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SAN_FLAGS}" )
endif()

option(WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
if(WARNINGS_AS_ERRORS)
  if(MSVC)
    add_compile_options(/WX)
  else()
    add_compile_options(-Werror)
  endif()
endif()

find_package(Eigen3 REQUIRED)

find_package(Threads)
if(USE_THREADS)
  find_package(TBB)
endif()

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

if(MSVC)
  add_compile_options(/wd4251)
endif()

if(WIN32)
    add_definitions( -DBOOST_ALL_NO_LIB
                     -DBOOST_ALL_DYN_LINK
                     -D_USE_MATH_DEFINES
                     -DNOMINMAX)
endif()

option(BUILD_EXAMPLES "Build SpinWaveGenie Examples" OFF)
option(BUILD_TESTING "Build unit tests" OFF)
option(ENABLE_DOXYGEN "Enables building of documentation with doxygen" OFF)
option(USE_THREADS "Use Threading Building Blocks to improve performance" ON)

if(BUILD_TESTING)
    enable_testing()
endif(BUILD_TESTING)

if(USE_THREADS AND TBB_FOUND )
    add_definitions(-DUSE_THREADS)
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
endif()

add_subdirectory(libSpinWaveGenie)

option( PYTHON_API "Build Python bindings for SpinWaveGenie" ON)
if(PYTHON_API)
  add_subdirectory(pythonapi)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

######################################################################
# configure everything for cpack
set ( CPACK_PACKAGE_NAME SpinWaveGenie )
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Library for simplifying linear spin wave calculations" )
set ( CPACK_PACKAGE_VERSION_MAJOR "${SPINWAVEGENIE_VERSION_MAJOR}" )
set ( CPACK_PACKAGE_VERSION_MINOR "${SPINWAVEGENIE_VERSION_MINOR}" )
set ( CPACK_PACKAGE_VERSION_PATCH "${SPINWAVEGENIE_VERSION_PATCH}" )
set ( CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )

set ( CPACK_PACKAGE_VENDOR "Steven Hahn" )
set ( CPACK_PACKAGE_URL "http://github.com/SpinWaveGenie/SpinWaveGenie" )
set ( CPACK_PACKAGE_CONTACT "hahnse@ornl.gov" )
set ( CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} )
set ( CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md )
set ( CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE )
set ( CPACK_RPM_PACKAGE_LICENSE "BSD 3-clause License" )

set ( CPACK_RPM_PACKAGE_RELEASE 1 )
set ( CPACK_RPM_PACKAGE_URL "${CPACK_PACKAGE_URL}" )
set ( CPACK_RPM_PACKAGE_GROUP "Development/Tools" )
set ( CPACK_RPM_PACKAGE_REQUIRES "eigen3-devel cmake >= 2.8.12 gcc-c++ doxygen boost-test boost-devel" )
if(USE_THREADS AND TBB_FOUND )
  set ( CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES} tbb-devel" )
  set ( CPACK_RPM_USE_THREADS "1" )
else()
  set ( CPACK_RPM_USE_THREADS "0" )
endif()

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/SpinWaveGenie.spec.in
                ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie.spec )

add_custom_command( OUTPUT SpinWaveGenie-${CPACK_PACKAGE_VERSION}.tar.gz
                    COMMAND git archive --format=tar --prefix=SpinWaveGenie-${CPACK_PACKAGE_VERSION}/ HEAD | gzip -c > ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie-${CPACK_PACKAGE_VERSION}.tar.gz
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    )
add_custom_target( source_tgz
                    DEPENDS SpinWaveGenie-${CPACK_PACKAGE_VERSION}.tar.gz )
