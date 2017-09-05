cmake_minimum_required(VERSION 2.8)
project(example_utopia_fe_dependent_project)

#path to the cmake_module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../cmake")

find_package(UtopiaFE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_FE_CXX_FLAGS}")
link_libraries(${UTOPIA_FE_LIBRARIES})
include_directories(${UTOPIA_FE_INCLUDES})
add_definitions(${UTOPIA_FE_DEFS})

add_executable(run_example ../utopia_fe_assembly_example.cpp)


