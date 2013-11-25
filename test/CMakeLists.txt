enable_testing()

include_directories(${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR}  ${spin_wave_genie_SOURCE_DIR}/src  .)

LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} cell)
add_executable(SublatticeTest SublatticeTest.cpp)
add_executable(CellTest CellTest.cpp) 
add_test(NAME SublatticeTest COMMAND SublatticeTest)
add_test(NAME CellTest COMMAND CellTest)


