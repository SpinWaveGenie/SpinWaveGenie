cmake_minimum_required(VERSION 3.0)

#add CFLAGS if not already present.
if(CMAKE_COMPILER_IS_GNUCXX AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+0x")
    message("compiler detected as: GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -msse4")
elseif(APPLE AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+11")
    message("Setting up Xcode for C++11 with libc++.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -msse4")
endif()

if(NOT SpinWaveGenie)
  find_package(SpinWaveGenie)
endif()
add_executable(YFeO3 TwoDimensionalCut.cpp CommonFunctions.cpp)
target_link_libraries( YFeO3 SpinWaveGenie)
target_include_directories( YFeO3 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file( COPY "YFeO3_2_K_m2.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "YFeO3_3_0_L.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "YFeO3_3_1_L.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
