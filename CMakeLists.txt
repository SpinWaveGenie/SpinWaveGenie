project(spin_wave_genie)

cmake_minimum_required (VERSION 2.8.12)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MACOSX_RPATH FALSE)
#set(CMAKE_SKIP_RPATH TRUE)

set(SPINWAVEGENIE_VERSION_MAJOR "0")
set(SPINWAVEGENIE_VERSION_MINOR "1")
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
    include(Coveralls)
    coveralls_turn_on_coverage()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "compiler detected as: GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE -std=c++0x -Wall -Wextra")
elseif(APPLE)
    message(STATUS "Setting up Xcode for C++11 with libc++.")
    if(CMAKE_GENERATOR STREQUAL Xcode)
        set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++0x")
        set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -msse4 -Wall -Wextra -Wpedantic")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libstdc++ -msse4 -Wall -Wextra -Wpedantic")
endif(CMAKE_COMPILER_IS_GNUCXX)

option(WITH_SANITIZERS "Enable sanitizers" OFF)
if(WITH_SANITIZERS)
  message(STATUS "enabling address sanitizer")
  set(SAN_FLAGS "-fno-omit-frame-pointer -fno-common -g -fsanitize=address,undefined -fno-sanitize-recover")
  add_compile_options(-fno-omit-frame-pointer -fno-common -g -fsanitize=address,undefined -fno-sanitize-recover)
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${SAN_FLAGS}" )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SAN_FLAGS}" )
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SAN_FLAGS}" )
endif()

find_package(Eigen3 REQUIRED)
find_package(TBB)

string(FIND ${CMAKE_SYSTEM} "el6" isRHEL6)
if (isRHEL6)
    set(Boost_NO_BOOST_CMAKE ON)
endif(isRHEL6)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

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

INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX( atomic HAVE_ATOMIC)
if(HAVE_ATOMIC)
    add_definitions(-DHAVE_ATOMIC_H)
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

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/SpinWaveGenie.spec.in
                 ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie.spec )


add_custom_command( OUTPUT SpinWaveGenie-${CPACK_PACKAGE_VERSION}.tar.gz
                    COMMAND git archive --format=tar --prefix=SpinWaveGenie-${CPACK_PACKAGE_VERSION}/ HEAD | gzip -c > ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie-${CPACK_PACKAGE_VERSION}.tar.gz
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    )
add_custom_target( source_tgz
                    DEPENDS SpinWaveGenie-${CPACK_PACKAGE_VERSION}.tar.gz )
