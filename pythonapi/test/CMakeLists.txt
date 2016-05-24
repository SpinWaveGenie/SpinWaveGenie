add_test(NAME python_SublatticeTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/SublatticeTest.py )
add_test(NAME python_CellTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/CellTest.py )
add_test(NAME python_NeighborsTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/NeighborsTest.py )
add_test(NAME python_InteractionFactoryTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/InteractionFactoryTest.py )
add_test(NAME python_SpinWaveBuilderTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/SpinWaveBuilderTest.py )

set_tests_properties(python_SublatticeTest python_CellTest python_NeighborsTest python_InteractionFactoryTest python_SpinWaveBuilderTest PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/../")
