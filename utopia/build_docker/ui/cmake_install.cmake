# Install script for directory: /shared/utopia/ui

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
    "/shared/utopia/ui/./utopia_AppRunner.hpp"
    "/shared/utopia/ui/./utopia_AppWithInputRegistry.hpp"
    "/shared/utopia/ui/./utopia_Convertible.hpp"
    "/shared/utopia/ui/./utopia_Input.hpp"
    "/shared/utopia/ui/./utopia_InputParameters.hpp"
    "/shared/utopia/ui/./utopia_JSONInput.hpp"
    "/shared/utopia/ui/./utopia_ObjectFactory.hpp"
    "/shared/utopia/ui/./utopia_SymbolicFunction.hpp"
    "/shared/utopia/ui/./utopia_TextOutput.hpp"
    "/shared/utopia/ui/./utopia_XMLInput.hpp"
    "/shared/utopia/ui/./utopia_ui.hpp"
    )
endif()

