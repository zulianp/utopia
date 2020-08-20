# Install script for directory: /shared/utopia/differentiable

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
    "/shared/utopia/differentiable/./utopia_AutoDiff.hpp"
    "/shared/utopia/differentiable/./utopia_AutoDiffExpr.hpp"
    "/shared/utopia/differentiable/./utopia_AutoDiffExpr_Binary.hpp"
    "/shared/utopia/differentiable/./utopia_AutoDiffExpr_Diag.hpp"
    "/shared/utopia/differentiable/./utopia_AutoDiffExpr_Multiply.hpp"
    "/shared/utopia/differentiable/./utopia_AutoDiffExpr_Reduce.hpp"
    "/shared/utopia/differentiable/./utopia_AutoDiffExpr_Trace.hpp"
    "/shared/utopia/differentiable/./utopia_AutoDiffExpr_Transposed.hpp"
    "/shared/utopia/differentiable/./utopia_Differentiable.hpp"
    "/shared/utopia/differentiable/./utopia_Simplify.hpp"
    "/shared/utopia/differentiable/./utopia_Simplify_Binary.hpp"
    "/shared/utopia/differentiable/./utopia_Simplify_Multiply.hpp"
    )
endif()

