enable_testing()

include_directories(${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR}  ${spin_wave_genie_SOURCE_DIR}/src  .)

LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES})

file(GLOB CELL_SRCS ../src/Cell/*.cpp)
file(GLOB CONTAINERS_SRCS ../src/Containers/*.cpp)
file(GLOB FORMFACTOR_SRCS ../src/Genie/MagneticFormFactor.cpp)

file(GLOB CELL_HDRS ../src/Cell/*.h)
file(GLOB CONTAINERS_HDRS ../src/Containers/*.h)
file(GLOB FORMFACTOR_HDRS ../src/Genie/MagneticFormFactor.h)

add_executable(SublatticeTest SublatticeTest.cpp ${CELL_SRCS} ${CELL_HDRS} ${CONTAINERS_SRCS} ${CONTAINERS_HDRS}  )
add_executable(CellTest CellTest.cpp ${CELL_SRCS} ${CELL_HDRS} ${CONTAINERS_HDRS} ${CONTAINERS_SRCS} ) 
add_executable(FormFactorTest FormFactorTest.cpp ${FORMFACTOR_SRCS} ${FORMFACTOR_HDRS} )
add_executable(HKLDirectionsTest HKLDirectionsTest.cpp ${CONTAINERS_SRCS} ${CONTAINERS_HDRS})
add_test(NAME SublatticeTest COMMAND SublatticeTest)
add_test(NAME CellTest COMMAND CellTest)
add_test(NAME FormFactorTest COMMAND FormFactorTest)
add_test(NAME HKLDirectionsTest COMMAND HKLDirectionsTest)


