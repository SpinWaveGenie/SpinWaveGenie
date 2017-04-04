cmake_minimum_required(VERSION 3.0)

if(NOT TARGET SpinWaveGenie)
  find_package(SpinWaveGenie)
endif()

add_executable(AFMDispersion Dispersion.cpp)
target_link_libraries( AFMDispersion SpinWaveGenie)
add_executable(AFMcut TwoDimensionalCut.cpp)
target_link_libraries( AFMcut SpinWaveGenie)
add_executable(AFMPowderAverage PowderAverage.cpp)
target_link_libraries( AFMPowderAverage SpinWaveGenie)

file( COPY "AFMChain.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "AFMavg.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

