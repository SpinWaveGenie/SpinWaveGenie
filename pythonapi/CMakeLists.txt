find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

if(NOT PYTHON_SITE_PACKAGES_DIR)
  execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(True))"
                    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  set( PYTHON_SITE_PACKAGES_DIR ${PYTHON_SITE_PACKAGES} CACHE FILEPATH "site-packages directory for python bindings")
endif()

message(STATUS "PYTHON_SITE_PACKAGES = ${PYTHON_SITE_PACKAGES_DIR}" )

pybind11_add_module( python_SpinWaveGenie src/SpinWaveGenie.cpp )
target_link_libraries( python_SpinWaveGenie PRIVATE SpinWaveGenie )

file(MAKE_DIRECTORY SpinWaveGenie)
set_target_properties ( python_SpinWaveGenie PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie )
set_target_properties ( python_SpinWaveGenie PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} )
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie)

if(CMAKE_GENERATOR MATCHES "Visual Studio")
  # Horrible hack so that the python library is found.
  set_target_properties ( python_SpinWaveGenie PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie
                                                          RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie
                                                          RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie
                                                          RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie )

elseif(CMAKE_GENERATOR MATCHES "Xcode")
  # Horrible hack so that the python library is found.
  set_target_properties ( python_SpinWaveGenie PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie
                                                          LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie
                                                          LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie
                                                          LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie )
endif()

if(WIN32)
  target_link_libraries( python_SpinWaveGenie SpinWaveGenie ${PYTHON_LIBRARY})
  set_target_properties( python_SpinWaveGenie PROPERTIES SUFFIX ".pyd")
else()
  set_target_properties( python_SpinWaveGenie PROPERTIES PREFIX "")
  set_target_properties( python_SpinWaveGenie PROPERTIES SUFFIX ".so")
endif()

if(APPLE)
  set_target_properties(python_SpinWaveGenie PROPERTIES LINK_FLAGS "-undefined dynamic_lookup ")
endif()

add_subdirectory(test)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${PYTHON_SITE_PACKAGES_DIR}/SpinWaveGenie)
install(TARGETS python_SpinWaveGenie DESTINATION ${PYTHON_SITE_PACKAGES_DIR}/SpinWaveGenie)

