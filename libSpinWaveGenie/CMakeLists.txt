file(GLOB EXTERNAL_SRCS src/External/*.cpp src/External/*.c)
file(GLOB CELL_SRCS src/Cell/*.cpp)
file(GLOB CONTAINERS_SRCS src/Containers/*.cpp)
file(GLOB GENIE_SRCS src/Genie/*.cpp)
file(GLOB INTERACTIONS_SRCS src/Interactions/*.cpp)
file(GLOB PLOT_SRCS src/Plot/*.cpp)
SOURCE_GROUP("Source Files\\External" FILES ${EXTERNAL_SRCS})
SOURCE_GROUP("Source Files\\Cell" FILES ${CELL_SRCS})
SOURCE_GROUP("Source Files\\Containers" FILES ${CONTAINERS_SRCS})
SOURCE_GROUP("Source Files\\Genie" FILES ${GENIE_SRCS})
SOURCE_GROUP("Source Files\\Interactions" FILES ${INTERACTIONS_SRCS})
SOURCE_GROUP("Source Files\\Plot" FILES ${PLOT_SRCS})

file(GLOB EXTERNAL_HDRS src/External/*.h)
file(GLOB CELL_HDRS include/SpinWaveGenie/Cell/*.h)
file(GLOB CONTAINERS_HDRS include/SpinWaveGenie/Containers/*.h)
file(GLOB GENIE_HDRS include/SpinWaveGenie/Genie/*.h)
file(GLOB INTERACTIONS_HDRS include/SpinWaveGenie/Interactions/*.h)
file(GLOB PLOT_HDRS include/SpinWaveGenie/Plot/*.h)
SOURCE_GROUP("Header Files\\External" FILES ${EXTERNAL_HDRS})
SOURCE_GROUP("Header Files\\Cell" FILES ${CELL_HDRS})
SOURCE_GROUP("Header Files\\Containers" FILES ${CONTAINERS_HDRS})
SOURCE_GROUP("Header Files\\Genie" FILES ${GENIE_HDRS})
SOURCE_GROUP("Header Files\\Interactions" FILES ${INTERACTIONS_HDRS})
SOURCE_GROUP("Header Files\\Plot" FILES ${PLOT_HDRS})
set(LIBRARY_SRCS ${EXTERNAL_SRCS} ${EXTERN_SRCS} ${CELL_SRCS} ${CONTAINERS_SRCS} ${GENIE_SRCS} ${INTERACTIONS_SRCS} ${PLOT_SRCS})
set(LIBRARY_HDRS ${EXTERNAL_HDRS} ${EXTERN_HDRS} ${CELL_HDRS} ${CONTAINERS_HDRS} ${GENIE_HDRS} ${INTERACTIONS_HDRS} ${PLOT_HDRS})

add_library(SpinWaveGenie SHARED ${LIBRARY_SRCS} ${LIBRARY_HDRS})

#We need need the header files both to build SpinWaveGenie and for the interface.
target_include_directories( SpinWaveGenie PRIVATE include )               
set_property( TARGET SpinWaveGenie PROPERTY
              INTERFACE_INCLUDE_DIRECTORIES
              "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
              "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>"
)
target_include_directories( SpinWaveGenie SYSTEM INTERFACE PUBLIC ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
target_include_directories( SpinWaveGenie PRIVATE src)
target_include_directories( SpinWaveGenie SYSTEM PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries( SpinWaveGenie PRIVATE ${TBB_LIBRARIES})

if(BUILD_TESTING)
    add_subdirectory(test)
endif (BUILD_TESTING)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS SpinWaveGenie EXPORT SpinWaveGenieTargets 
                              ARCHIVE DESTINATION lib 
                              LIBRARY DESTINATION lib)

#create package configuration files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfigVersion.cmake"
  VERSION ${SPINWAVEGENIE_VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT SpinWaveGenieTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieTargets.cmake"
  NAMESPACE SpinWaveGenie::
)

configure_file(SpinWaveGenieConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfig.cmake"
  COPYONLY
)

#install package configuration files
set(ConfigPackageLocation lib/cmake/SpinWaveGenie)
install(EXPORT SpinWaveGenieTargets
  FILE
    SpinWaveGenieTargets.cmake
  NAMESPACE
    Upstream::
  DESTINATION
    ${ConfigPackageLocation}
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)

string(MD5 md5_PackageLocation "${CMAKE_INSTALL_PREFIX}/${ConfigPackageLocation}/SpinWaveGenieConfig.cmake")
file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/${md5_PackageLocation} "${CMAKE_INSTALL_PREFIX}/${ConfigPackageLocation}/SpinWaveGenieConfig.cmake")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${md5_PackageLocation}" DESTINATION "$ENV{HOME}/.cmake/packages/SpinWaveGenie/")

