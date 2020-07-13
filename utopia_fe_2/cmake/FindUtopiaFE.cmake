# -*- mode: cmake -*-
#
# UTOPIA_FE Find Module
#
# Usage:
#    Control the search through UTOPIA_FE_DIR.
#
#    Following variables are set:
#    UTOPIA_FE_FOUND            (BOOL)               Flag indicating if UTOPIA was found
#    UTOPIA_FE_INCLUDES         (LIST of PATH)       Path to the UTOPIA include file
#    UTOPIA_FE_DEFS             (LIST)               List with all defs for the compiler
#    UTOPIA_FE_LIBRARIES        (LIST)               List of all required UTOPIA libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Search for include files
find_path(UTOPIA_FE_INSTALLATION_PATH
  NAMES "include/utopia_fe.hpp"
  HINTS
      $ENV{UTOPIA_FE_DIR}
      ${UTOPIA_FE_DIR}
      $ENV{UTOPIA_DIR}
      $UTOPIA_DIR
  NO_DEFAULT_PATH)


find_path(UTOPIA_FE_CONFIG_FILE_PATH
  NAMES "utopia_fe_config.cmake"
  HINTS
    $ENV{UTOPIA_FE_DIR}
    ${UTOPIA_FE_DIR}
    ${UTOPIA_FE_INCLUDES}/..
    $ENV{UTOPIA_DIR}
    ${UTOPIA_DIR}
    PATH_SUFFIXES config
  NO_DEFAULT_PATH)

if(UTOPIA_FE_CONFIG_FILE_PATH)
  include(${UTOPIA_FE_CONFIG_FILE_PATH}/utopia_fe_config.cmake)
endif()

# Send useful message if everything is found
find_package_handle_standard_args(
  UTOPIA DEFAULT_MSG
  UTOPIA_FE_LIBRARIES
  UTOPIA_FE_INCLUDES)

# find_package_handle_standard_args should set UTOPIA_FE_FOUND but it does not!
if (UTOPIA_FE_LIBRARIES AND UTOPIA_FE_INCLUDES AND UTOPIA_FE_CONFIG_FILE_PATH)
  set(UTOPIA_FE_FOUND TRUE)
else()
  set(UTOPIA_FE_FOUND FALSE)
  message(WARNING "utopia not found. You can set UTOPIA_FE_DIR variable to the path of the utopia installation. Either in your environment or through -DUTOPIA_FE_DIR=...")
  MESSAGE(STATUS "---------------------------------------------")
  MESSAGE(STATUS "utopia_fe include: ${UTOPIA_FE_INCLUDES}")
  MESSAGE(STATUS "utopia_fe lib:     ${UTOPIA_FE_LIBRARIES}")
  MESSAGE(STATUS "utopia_fe config:  ${UTOPIA_FE_CONFIG_FILE_PATH}")
  MESSAGE(STATUS "---------------------------------------------")
endif()


mark_as_advanced(UTOPIA_FE_INCLUDES UTOPIA_FE_LIBRARIES UTOPIA_FE_CONFIG_FILE_PATH)

