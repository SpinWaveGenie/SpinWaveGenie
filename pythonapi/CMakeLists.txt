#execute_process(COMMAND python-config --prefix OUTPUT_VARIABLE PYTHON_PREFIX)
#string(STRIP ${PYTHON_PREFIX} PYTHON_PREFIX)
#set( PYTHON_LIBRARY "${PYTHON_PREFIX}/lib/libpython2.7.dylib" CACHE FILEPATH "PYTHON_LIBRARY")
#set( PYTHON_INCLUDE_DIR "${PYTHON_PREFIX}/include/python2.7" CACHE PATH "PYTHON_INCLUDE_DIR")

find_package ( PythonLibs REQUIRED )

include(ExternalProject)
ExternalProject_Add(pybind11 GIT_REPOSITORY https://github.com/pybind/pybind11.git GIT_TAG master PREFIX ${CMAKE_CURRENT_BINARY_DIR}/pybind11 BUILD_COMMAND "" INSTALL_COMMAND "")
ExternalProject_Get_Property(pybind11 source_dir)

add_library( PySpinWaveGenie SHARED src/SpinWaveGenie.cpp)

target_include_directories( PySpinWaveGenie PRIVATE ${PYTHON_INCLUDE_DIR} ${source_dir}/include )
target_link_libraries ( PySpinWaveGenie SpinWaveGenie ${PYTHON_LIBRARY})
set_target_properties( PySpinWaveGenie PROPERTIES PREFIX "")
set_target_properties( PySpinWaveGenie PROPERTIES SUFFIX ".so")

add_subdirectory(test)
