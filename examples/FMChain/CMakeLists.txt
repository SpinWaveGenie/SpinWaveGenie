cmake_minimum_required(VERSION 3.0)

if(NOT TARGET SpinWaveGenie)
  find_package(SpinWaveGenie)
endif()
add_executable(FMDispersion Dispersion.cpp)
target_link_libraries( FMDispersion SpinWaveGenie)
add_executable(FMcut TwoDimensionalCut.cpp)
target_link_libraries( FMcut SpinWaveGenie)
add_executable(FMPowderAverage PowderAverage.cpp)
target_link_libraries( FMPowderAverage SpinWaveGenie)

file( COPY "FM.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "FMavg.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
