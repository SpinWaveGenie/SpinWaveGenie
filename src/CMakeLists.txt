include_directories(${EIGEN3_INCLUDE_DIR})
set(SRCS main.cpp SW_sublattice.cpp SW_Matrix.cpp)
add_executable(../spin_wave_genie ${SRCS})
