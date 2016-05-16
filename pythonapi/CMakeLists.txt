execute_process(COMMAND python-config --prefix OUTPUT_VARIABLE PYTHON_PREFIX)
string(STRIP ${PYTHON_PREFIX} PYTHON_PREFIX)
set( PYTHON_LIBRARY "${PYTHON_PREFIX}/lib/libpython2.7.dylib" CACHE FILEPATH "PYTHON_LIBRARY")
set( PYTHON_INCLUDE_DIR "${PYTHON_PREFIX}/include/python2.7" CACHE PATH "PYTHON_INCLUDE_DIR")

find_package ( PythonLibs REQUIRED )

add_library( PySpinWaveGenie SHARED SpinWaveGenie.cpp)

target_include_directories( PySpinWaveGenie PRIVATE ${PYTHON_INCLUDE_DIR} /usr/local/include )
target_link_libraries ( PySpinWaveGenie SpinWaveGenie ${PYTHON_LIBRARY})
set_target_properties( PySpinWaveGenie PROPERTIES PREFIX "")
set_target_properties( PySpinWaveGenie PROPERTIES SUFFIX ".so")

