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
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie)
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  add_custom_command ( TARGET python_SpinWaveGenie POST_BUILD
                       COMMAND ${CMAKE_COMMAND} ARGS -E copy
                       ${CMAKE_CFG_INTDIR}/python_SpinWaveGenie.pyd
                       SpinWaveGenie/python_SpinWaveGenie.pyd )
else ()
  set_target_properties ( python_SpinWaveGenie PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie )
endif ()

set_target_properties ( python_SpinWaveGenie PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} )

add_subdirectory(test)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${PYTHON_SITE_PACKAGES_DIR}/SpinWaveGenie)
install(TARGETS python_SpinWaveGenie DESTINATION ${PYTHON_SITE_PACKAGES_DIR}/SpinWaveGenie)

