cmake_minimum_required(VERSION 2.8)
project(example_utopia_fe_dependent_project)

#path to the cmake module FindUtopiaFE.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../cmake")

find_package(UtopiaFE)

if(UTOPIA_FE_FOUND)
    add_executable(run_example ../utopia_fe_assembly_example.cpp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_FE_CXX_FLAGS}")
    target_link_libraries(run_example UtopaFE::utopia_fe)
else()
    message(FATAL_ERROR "[Problem] Cannot compile project without utopia_fe.\n[Solution] Set the variable UTOPIA_FE_DIR to the installation location of utopia_fe (i.e., the install prefix)")
endif()
