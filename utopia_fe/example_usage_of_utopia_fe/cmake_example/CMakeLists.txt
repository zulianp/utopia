cmake_minimum_required(VERSION 3.13)
project(example_utopia_fe_dependent_project)

find_package(UtopiaFE REQUIRED)
add_executable(run_example ../utopia_fe_assembly_example.cpp)
target_link_libraries(run_example UtopiaFE::utopia_fe)
