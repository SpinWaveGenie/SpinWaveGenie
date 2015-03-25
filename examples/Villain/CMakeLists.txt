cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#add CFLAGS if not already present.
if(CMAKE_COMPILER_IS_GNUCXX AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+0x")
    message("compiler detected as: GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -msse4")
elseif(APPLE AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+11")
    message("Setting up Xcode for C++11 with libc++.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -msse4")
endif()

find_package(NLOPT)
find_package(SpinWaveGenie)
add_executable(VillainSC2Dispersion SC2Dispersion.cpp)
target_link_libraries( VillainSC2Dispersion SpinWaveGenie)
add_executable(VillainSC1Dispersion SC1Dispersion.cpp)
target_link_libraries( VillainSC1Dispersion SpinWaveGenie ${NLOPT_LIBRARIES})
target_include_directories(VillainSC1Dispersion SYSTEM PRIVATE ${NLOPT_INCLUDE_DIRECTORY})
add_executable(VillainFMDispersion FMDispersion.cpp)
target_link_libraries( VillainFMDispersion SpinWaveGenie)
