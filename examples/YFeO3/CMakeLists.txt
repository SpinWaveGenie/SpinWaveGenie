cmake_minimum_required(VERSION 3.0)

if(NOT TARGET SpinWaveGenie)
  find_package(SpinWaveGenie)
endif()

add_executable(YFeO3 TwoDimensionalCut.cpp CommonFunctions.cpp)
target_link_libraries( YFeO3 SpinWaveGenie)
target_include_directories( YFeO3 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file( COPY "YFeO3_2_K_m2.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "YFeO3_2_K_m3.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "YFeO3_3_0_L.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "YFeO3_3_1_L.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
