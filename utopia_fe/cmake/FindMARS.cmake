# -*- mode: cmake -*-
#
# MARS Find Module
#
# Usage:
#    Control the search through MARS_DIR.
#
#    Following variables are set:
#    MARS_FOUND            (BOOL)               Flag indicating if MARS was found
#    MARS_INCLUDES         (LIST of PATH)       Path to the MARS include file
#    MARS_DEFS             (LIST)               List with all defs for the compiler
#    MARS_LIBRARIES        (LIST)               List of all required MARS libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Search for include files
find_path(MARS_INCLUDES
  NAMES "mars.hpp"
  HINTS 
      ${MARS_DIR}
      $ENV{MARS_DIR}
  PATH_SUFFIXES include  
  NO_DEFAULT_PATH)

# Search for libraries 
find_library(MARS_LIBRARIES
 NAMES mars
 HINTS 
     ${MARS_DIR}
     ${MARS_INCLUDES}/..
     $ENV{MARS_DIR}
     PATH_SUFFIXES lib 
 NO_DEFAULT_PATH)

# find_path(MARS_CONFIG_FILE_PATH 
#   NAMES "mars_config.cmake"
#   HINTS
#     $ENV{MARS_DIR}
#     ${MARS_DIR}
#     ${MARS_INCLUDES}/..
#     $ENV{MOONOLITH_ROOT}/utopia/build
#     PATH_SUFFIXES config 
#   NO_DEFAULT_PATH)



# if(MARS_CONFIG_FILE_PATH)
#   include(${MARS_CONFIG_FILE_PATH}/mars_config.cmake)
# endif()

# Send useful message if everything is found
find_package_handle_standard_args(
  MARS DEFAULT_MSG
  MARS_LIBRARIES
  MARS_INCLUDES)

# find_package_handle_standard_args should set MARS_FOUND but it does not!
if (MARS_LIBRARIES AND MARS_INCLUDES)# AND MARS_CONFIG_FILE_PATH)
  set(MARS_FOUND TRUE)
else()
  set(MARS_FOUND FALSE)
  message(WARNING "mars not found. You can set MARS_DIR variable to the path of the utopia installation. Either in your environment or through -DMARS_DIR=...")
  MESSAGE(STATUS "---------------------------------------------")
  MESSAGE(STATUS "mars include: ${MARS_INCLUDES}")
  MESSAGE(STATUS "mars lib:     ${MARS_LIBRARIES}")
  MESSAGE(STATUS "mars config:  ${MARS_CONFIG_FILE_PATH}")
  MESSAGE(STATUS "---------------------------------------------")
endif()

mark_as_advanced(MARS_INCLUDES MARS_LIBRARIES MARS_CONFIG_FILE_PATH)

