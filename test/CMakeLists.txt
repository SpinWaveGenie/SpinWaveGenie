enable_testing()

include_directories(${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR}  ${spin_wave_genie_SOURCE_DIR}/src  .)

LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} cell)
add_executable(test_sublattice test_sublattice.cpp)
 
add_test(NAME test_sublattice COMMAND test_sublattice)



