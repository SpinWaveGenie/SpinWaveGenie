include_directories(${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR} SpinWaveGenie ${spin_wave_genie_SOURCE_DIR}/libSpinWaveGenie/include)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} SpinWaveGenie)

add_executable(AFMDispersion AFMChain/Dispersion.cpp)
add_executable(AFMcut AFMChain/TwoDimensionalCut.cpp)

add_executable(FMPowderAverage FMChain/PowderAverage.cpp)
add_executable(FMDispersion FMChain/Dispersion.cpp)
add_executable(FMcut FMChain/TwoDimensionalCut.cpp)
