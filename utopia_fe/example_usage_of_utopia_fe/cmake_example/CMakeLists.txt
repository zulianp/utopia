cmake_minimum_required(VERSION 3.13)
project(example_utopia_fe_dependent_project)

include(CMakeFindDependencyMacro)

find_package(UtopiaFE REQUIRED)

if(UtopiaFE_FOUND)
  message(STATUS "UtopiaFE Found!")
else()
  message(STATUS "UtopiaFE Not Found!")
endif()

add_executable(run_example ../utopia_fe_assembly_example.cpp)
target_link_libraries(run_example UtopiaFE::utopia_fe)

enable_testing()

add_test(NAME UTOPIA_FE.utopia_fe_post_install_example_1 COMMAND run_example)