# -*- mode: cmake -*-
#
# UTOPIA Find Module
#
# Usage:
#    Control the search through UTOPIA_DIR.
#
#    Following variables are set:
#    UTOPIA_FOUND            (BOOL)               Flag indicating if UTOPIA was found
#    UTOPIA_INCLUDES         (LIST of PATH)       Path to the UTOPIA include file
#    UTOPIA_DEFS             (LIST)               List with all defs for the compiler
#    UTOPIA_LIBRARIES        (LIST)               List of all required UTOPIA libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Search for include files
find_path(UTOPIA_INCLUDES
  NAMES "utopia.hpp"
  HINTS
      ${UTOPIA_DIR}
      $ENV{UTOPIA_DIR}
      $ENV{MOONOLITH_ROOT}/utopia/build
  PATH_SUFFIXES include
  NO_DEFAULT_PATH)

# Search for libraries
find_library(UTOPIA_LIBRARIES
 NAMES utopia
 HINTS
     ${UTOPIA_DIR}
     ${UTOPIA_INCLUDES}/..
     $ENV{UTOPIA_DIR}

     PATH_SUFFIXES lib
 NO_DEFAULT_PATH)

find_path(UTOPIA_CONFIG_FILE_PATH
  NAMES "utopia-config.cmake"
  HINTS
    ${UTOPIA_DIR}
    $ENV{UTOPIA_DIR}
    ${UTOPIA_INCLUDES}/..
    $ENV{MOONOLITH_ROOT}/utopia/build
    PATH_SUFFIXES config
  NO_DEFAULT_PATH)



if(UTOPIA_CONFIG_FILE_PATH)
  include(${UTOPIA_CONFIG_FILE_PATH}/utopia-config.cmake)
endif()

# Send useful message if everything is found
find_package_handle_standard_args(
  UTOPIA DEFAULT_MSG
  UTOPIA_LIBRARIES
  UTOPIA_INCLUDES)

# find_package_handle_standard_args should set UTOPIA_FOUND but it does not!
if (UTOPIA_LIBRARIES AND UTOPIA_INCLUDES AND UTOPIA_CONFIG_FILE_PATH)
  set(UTOPIA_FOUND TRUE)
else()
  set(UTOPIA_FOUND FALSE)
  message(WARNING "utopia not found. You can set UTOPIA_DIR variable to the path of the utopia installation. Either in your environment or through -DUTOPIA_DIR=...")
  MESSAGE(STATUS "---------------------------------------------")
  MESSAGE(STATUS "moonolith root: $ENV{MOONOLITH_ROOT}")
  MESSAGE(STATUS "utopia include: ${UTOPIA_INCLUDES}")
  MESSAGE(STATUS "utopia lib:     ${UTOPIA_LIBRARIES}")
  MESSAGE(STATUS "utopia config:  ${UTOPIA_CONFIG_FILE_PATH}")
  MESSAGE(STATUS "---------------------------------------------")
endif()



mark_as_advanced(UTOPIA_INCLUDES UTOPIA_LIBRARIES UTOPIA_CONFIG_FILE_PATH)

