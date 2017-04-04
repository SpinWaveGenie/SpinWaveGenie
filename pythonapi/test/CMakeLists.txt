add_test(NAME python_SublatticeTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/SublatticeTest.py )
add_test(NAME python_CellTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/CellTest.py )
add_test(NAME python_NeighborsTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/NeighborsTest.py )
add_test(NAME python_InteractionFactoryTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/InteractionFactoryTest.py )
add_test(NAME python_FMDispersionTest COMMAND ${PYTHON_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/FMDispersionTest.py )

set(PYTHON_TEST_NAMES python_SublatticeTest
                      python_CellTest
                      python_NeighborsTest
                      python_InteractionFactoryTest
                      python_FMDispersionTest )

if (CMAKE_GENERATOR MATCHES "Visual Studio" OR CMAKE_GENERATOR MATCHES "Xcode")
  set_tests_properties( ${PYTHON_TEST_NAMES}
                        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/../$<CONFIGURATION>" )
else ()
  set_tests_properties( ${PYTHON_TEST_NAMES} 
                        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/../")
endif ()
