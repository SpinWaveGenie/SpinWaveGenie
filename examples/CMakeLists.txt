include_directories( ${EIGEN3_INCLUDE_DIR} ${spin_wave_genie_SOURCE_DIR}/libSpinWaveGenie/include ${NLOPT_INCLUDE_DIRECTORY})
LINK_DIRECTORIES()
LINK_LIBRARIES(SpinWaveGenie ${NLOPT_LIBRARIES})

add_executable(AFMDispersion AFMChain/Dispersion.cpp)
add_executable(AFMcut AFMChain/TwoDimensionalCut.cpp)
add_executable(AFMPowderAverage AFMChain/PowderAverage.cpp)

add_executable(FMPowderAverage FMChain/PowderAverage.cpp)
add_executable(FMDispersion FMChain/Dispersion.cpp)
add_executable(FMcut FMChain/TwoDimensionalCut.cpp)

add_executable(VillainSC2Dispersion Villain/SC2Dispersion.cpp)
add_executable(VillainSC1Dispersion Villain/SC1Dispersion.cpp)
add_executable(VillainFMDispersion  Villain/FMDispersion.cpp)
