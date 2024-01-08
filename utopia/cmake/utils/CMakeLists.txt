cmake_minimum_required(VERSION 3.5)
project(utopia_gen_makefile_config)

include(CMakeFindDependencyMacro)

find_package(Utopia REQUIRED)

include(../GenerateMakeConfig.cmake)

generate_make_config(utopia Utopia::utopia)
