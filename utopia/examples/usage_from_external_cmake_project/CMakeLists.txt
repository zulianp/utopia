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

enable_testing()
# add_test(NAME UTOPIA.utopia_post_install_example_1 COMMAND utopia_post_install_example_1)
