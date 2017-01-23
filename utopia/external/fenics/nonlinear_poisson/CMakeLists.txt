#  cmake .. -DUTOPIA_DIR=/Users/alenakopanicakova/MOOSE/moonolith/utopia_installation
#
#
# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME demo_nonlinear-poisson)
project(${PROJECT_NAME})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../../../cmakemodules")
message(STATUS "${CMAKE_MODULE_PATH}")

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in
# DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN)

# Need to get VTK config because VTK uses advanced VTK features which
# mean it's not enough to just link to the DOLFIN target. See
# http://www.vtk.org/pipermail/vtk-developers/2013-October/014402.html
find_package(VTK HINTS ${VTK_DIR} $ENV{VTK_DIR} NO_MODULE QUIET)

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()



find_package(Utopia)
if(NOT UTOPIA_FOUND)
	message(ERROR "needs utopia")
endif()


# Add utopia library 
include("$ENV{UTOPIA_DIR}/config/utopia_config.cmake")

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})
add_definitions(${UTOPIA_DEFS})
# add_definitions(${CMAKE_C_COMPILER})
# add_definitions(${CMAKE_CXX_COMPILER})

# Compiler flags
set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
list(APPEND ${CMAKE_CXX_FLAGS})



# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS})
include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})



# Utopia dirs 
include_directories(${UTOPIA_INCLUDES})

# Executable
add_executable(${PROJECT_NAME} main.cpp)

# Target libraries
target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES} ${UTOPIA_LIBRARIES})
# target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES})

message(STATUS "${UTOPIA_LIBRARIES}")
