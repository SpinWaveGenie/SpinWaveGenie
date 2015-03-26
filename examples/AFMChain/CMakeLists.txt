cmake_minimum_required(VERSION 3.0)

#add CFLAGS if not already present.
if(CMAKE_COMPILER_IS_GNUCXX AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+0x")
    message(STATUS "compiler detected as: GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -msse4")
elseif(APPLE AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+11")
    message(STATUS "Setting up Xcode for C++11 with libc++.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -msse4")
endif()

find_package(SpinWaveGenie)
add_executable(AFMDispersion Dispersion.cpp)
target_link_libraries( AFMDispersion SpinWaveGenie)
add_executable(AFMcut TwoDimensionalCut.cpp)
target_link_libraries( AFMcut SpinWaveGenie)
add_executable(AFMPowderAverage PowderAverage.cpp)
target_link_libraries( AFMPowderAverage SpinWaveGenie)

file( COPY "AFMChain.py"
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
