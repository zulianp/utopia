cmake_minimum_required(VERSION 3.5)
project(utopia_fe_gen_makefile_config)

include(CMakeFindDependencyMacro)

find_package(UtopiaFE REQUIRED)

include(../GenerateMakeConfig.cmake)

generate_make_config(utopia_fe UtopiaFE::utopia_fe)
