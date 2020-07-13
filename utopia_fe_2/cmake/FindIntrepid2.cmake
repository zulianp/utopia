#
# Try to find the Intrepid2 library
#
# This module exports:
#
#   INTREPID2_DIR
#   INTREPID2_INCLUDE_DIRS
#   INTREPID2_LIBRARIES
#   INTREPID2_VERSION
#   INTREPID2_VERSION_MAJOR
#   INTREPID2_VERSION_MINOR
#   INTREPID2_VERSION_SUBMINOR
#   INTREPID2_WITH_MANDATORY_CXX11
#   INTREPID2_WITH_MPI
#   INTREPID2_SUPPORTS_CPP11
#   INTREPID2_HAS_C99_TR1_WORKAROUND
#

#
# Include the trilinos package configuration:
#
FIND_PACKAGE(INTREPID2_CONFIG
  CONFIG
  NAMES Intrepid2 Intrepid2Config.cmake Intrepid2Config
  HINTS
    ${TRILINOS_DIR}/lib/cmake/Intrepid2
    ${INTREPID2_DIR}
    ${INTREPID2_DIR}
    $ENV{INTREPID2_DIR}
    $ENV{TRILINOS_DIR}
  PATH_SUFFIXES
    lib64/cmake/Intrepid2
    lib/cmake/Intrepid2
    lib${LIB_SUFFIX}/cmake/Intrepid2
  )

IF(DEFINED Intrepid2_VERSION)
  #
  # Extract version numbers:
  #
  SET(INTREPID2_VERSION "${Intrepid2_VERSION}")

  STRING(REGEX REPLACE
    "^([0-9]+).*$" "\\1"
    INTREPID2_VERSION_MAJOR "${Intrepid2_VERSION}")

  STRING(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*$" "\\1"
    INTREPID2_VERSION_MINOR "${Intrepid2_VERSION}")

  # If there is no subminor number,
  # INTREPID2_VERSION_SUBMINOR is set to an empty string.
  # If that is the case, set the subminor number to zero
  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.?(([0-9]+)?).*$" "\\1"
    INTREPID2_VERSION_SUBMINOR "${Intrepid2_VERSION}")
  IF("${INTREPID2_VERSION_SUBMINOR}" STREQUAL "")
    SET(INTREPID2_VERSION_SUBMINOR "0")
  ENDIF()
ENDIF()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Intrepid2 DEFAULT_MSG  Intrepid2_INCLUDE_DIRS Intrepid2_LIBRARIES Intrepid2_TPL_LIBRARIES)

