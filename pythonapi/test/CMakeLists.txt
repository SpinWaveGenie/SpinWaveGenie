find_package(PythonInterp REQUIRED)

add_test(NAME PySublatticeTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/SublatticeTest.py )
add_test(NAME PyCellTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/CellTest.py )
add_test(NAME PyNeighborsTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/NeighborsTest.py )

set_tests_properties(PySublatticeTest PyCellTest PyNeighborsTest PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/../")

