include_directories(${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR}  Interactions .)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} /usr/local/Cellar/nlopt/2.4/lib)
LINK_LIBRARIES(${Boost_LIBRARIES} nlopt)

#add_subdirectory(Interactions)

file(GLOB SRCS *.cpp *.h)
file(GLOB EXTERN_SRCS extern/*.cpp extern/*.c extern/*.h)
file(GLOB CELL_SRCS Cell/*.cpp Cell/*.h)
file(GLOB INTERACTIONS_SRCS Interactions/*.cpp Interactions/*.h)

set(ALL_SRCS ${CELL_SRCS} ${EXTERN_SRCS} ${SRCS} ${INTERACTIONS_SRCS}) 

add_executable(spin_wave_genie ${ALL_SRCS})
add_library(cell SHARED ${CELL_SRCS})
