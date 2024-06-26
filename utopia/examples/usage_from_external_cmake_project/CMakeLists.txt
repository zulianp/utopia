cmake_minimum_required(VERSION 3.5...3.19)
project(utopia_usage LANGUAGES "CXX" "C")

include(CMakeFindDependencyMacro)
find_package(Utopia REQUIRED)

if(Utopia_FOUND)
  message(STATUS "Utopia Found!")
else()
  message(STATUS "Utopia Not Found!")
endif()

add_executable(utopia_post_install_example_1 ${CMAKE_SOURCE_DIR}/../example_hello_world.cpp)
target_link_libraries(utopia_post_install_example_1 Utopia::utopia)

add_executable(utopia_post_install_example_2 ${CMAKE_SOURCE_DIR}/../example_communicator.cpp)
target_link_libraries(utopia_post_install_example_2 Utopia::utopia)

add_executable(utopia_post_install_example_3 ${CMAKE_SOURCE_DIR}/../example_comparison_linear_solvers.cpp)
target_link_libraries(utopia_post_install_example_3 Utopia::utopia)

add_executable(utopia_post_install_example_4 ${CMAKE_SOURCE_DIR}/../example_linear_regression.cpp)
target_link_libraries(utopia_post_install_example_4 Utopia::utopia)

add_executable(utopia_post_install_example_5 ${CMAKE_SOURCE_DIR}/../example_loops.cpp)
target_link_libraries(utopia_post_install_example_5 Utopia::utopia)

add_executable(utopia_post_install_example_6 ${CMAKE_SOURCE_DIR}/../example_nonlinear_solve.cpp)
target_link_libraries(utopia_post_install_example_6 Utopia::utopia)

add_executable(utopia_post_install_example_7 ${CMAKE_SOURCE_DIR}/../example_solve_PDE_1D.cpp)
target_link_libraries(utopia_post_install_example_7 Utopia::utopia)

add_executable(utopia_post_install_example_8 ${CMAKE_SOURCE_DIR}/../example_sum_of_two_vectors.cpp)
target_link_libraries(utopia_post_install_example_8 Utopia::utopia)

add_executable(utopia_post_install_example_9 ${CMAKE_SOURCE_DIR}/../example_trivial_linear_system.cpp)
target_link_libraries(utopia_post_install_example_9 Utopia::utopia)


enable_testing()
add_test(NAME UTOPIA.utopia_post_install_example_1 COMMAND utopia_post_install_example_1)
add_test(NAME UTOPIA.utopia_post_install_example_2 COMMAND utopia_post_install_example_2)
add_test(NAME UTOPIA.utopia_post_install_example_3 COMMAND utopia_post_install_example_3)
add_test(NAME UTOPIA.utopia_post_install_example_5 COMMAND utopia_post_install_example_5)
add_test(NAME UTOPIA.utopia_post_install_example_6 COMMAND utopia_post_install_example_6)
add_test(NAME UTOPIA.utopia_post_install_example_7 COMMAND utopia_post_install_example_7)
add_test(NAME UTOPIA.utopia_post_install_example_8 COMMAND utopia_post_install_example_8)
add_test(NAME UTOPIA.utopia_post_install_example_9 COMMAND utopia_post_install_example_9)
