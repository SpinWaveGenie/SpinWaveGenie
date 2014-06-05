enable_testing()

include_directories(${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR}  ${spin_wave_genie_SOURCE_DIR}/src  .)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} SpinWaveGenie)

add_executable(SublatticeTest SublatticeTest.cpp)
add_executable(CellTest CellTest.cpp ) 
add_executable(FormFactorTest FormFactorTest.cpp )
add_executable(HKLDirectionsTest HKLDirectionsTest.cpp )
add_executable(EnergiesTest EnergiesTest.cpp )
add_executable(ResultsTest ResultsTest.cpp )

add_test(NAME SublatticeTest COMMAND SublatticeTest)
add_test(NAME CellTest COMMAND CellTest)
add_test(NAME FormFactorTest COMMAND FormFactorTest)
add_test(NAME HKLDirectionsTest COMMAND HKLDirectionsTest)
add_test(NAME EnergiesTest COMMAND EnergiesTest)
add_test(NAME ResultsTest COMMAND ResultsTest)




