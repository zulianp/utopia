# This file is automatically generated by running
#
#  cmake .. -DUTOPIA_DIR=/Users/alenakopanicakova/MOOSE/moonolith/utopia_installation
#     cmake/scripts/generate-cmakefiles
#
# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME nonlinear_elasticity)
project(${PROJECT_NAME})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../../../cmakemodules")
message(STATUS "${CMAKE_MODULE_PATH}")

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in
# DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)
find_package(Utopia)
if(NOT UTOPIA_FOUND)
  message(ERROR "needs utopia")
endif()


if (EXISTS ${DOLFIN_USE_FILE})
  include(${DOLFIN_USE_FILE})

  # Default build type (can be overridden by user)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
  endif()
else()
  # Compiler definitions
  include("../../../../utopia_installation/config/utopia-config.cmake")
  add_definitions(${DOLFIN_CXX_DEFINITIONS})
  add_definitions(${UTOPIA_DEFS})

  # Compiler flags
  set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  list(APPEND ${CMAKE_CXX_FLAGS})

  # Include directories
  include_directories(${DOLFIN_INCLUDE_DIRS})
  include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
endif()


include_directories(${UTOPIA_INCLUDES})
# Executable
add_executable(${PROJECT_NAME} main.cpp)

# Target libraries
target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES} ${UTOPIA_LIBRARIES})

message(STATUS "${UTOPIA_LIBRARIES}")
