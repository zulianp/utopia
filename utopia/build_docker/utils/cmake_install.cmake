# Install script for directory: /shared/utopia/utils

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/shared/utopia/utils/./utopia_AuthoredWork.hpp"
    "/shared/utopia/utils/./utopia_CSV.hpp"
    "/shared/utopia/utils/./utopia_Chrono.hpp"
    "/shared/utopia/utils/./utopia_CiteUtopia.hpp"
    "/shared/utopia/utils/./utopia_ConditionalType.hpp"
    "/shared/utopia/utils/./utopia_Describable.hpp"
    "/shared/utopia/utils/./utopia_MPI.hpp"
    "/shared/utopia/utils/./utopia_Path.hpp"
    "/shared/utopia/utils/./utopia_Utils.hpp"
    "/shared/utopia/utils/./utopia_make_unique.hpp"
    "/shared/utopia/utils/./utopia_std_function.hpp"
    "/shared/utopia/utils/action/utopia_ActionMacros.hpp"
    "/shared/utopia/utils/action/utopia_ActionRegistry.hpp"
    "/shared/utopia/utils/action/utopia_NaryActionRegistry.hpp"
    "/shared/utopia/utils/action/utopia_NaryActionRegistry_impl.hpp"
    "/shared/utopia/utils/unit_testing/utopia_TestMacros.hpp"
    "/shared/utopia/utils/unit_testing/utopia_TestRegistry.hpp"
    "/shared/utopia/utils/unit_testing/utopia_TestRunner.hpp"
    "/shared/utopia/utils/unit_testing/utopia_Testing.hpp"
    "/shared/utopia/utils/app_management/utopia_AppMacros.hpp"
    "/shared/utopia/utils/app_management/utopia_AppRegistry.hpp"
    )
endif()

