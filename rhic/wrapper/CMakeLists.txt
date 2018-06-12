cmake_minimum_required (VERSION 2.6)
project (testHydroWrapper)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

set (CMAKE_CXX_FLAGS "-g -O3 -Wall")

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_executable(testHydroWrapper testHydroWrapper.cpp)

target_link_libraries(testHydroWrapper ${LIBS})
