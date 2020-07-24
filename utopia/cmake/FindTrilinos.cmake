#
# Try to find the Trilinos library
#
# This module exports:
#
#   TRILINOS_DIR
#   TRILINOS_INCLUDE_DIRS
#   TRILINOS_LIBRARIES
#   TRILINOS_VERSION
#   TRILINOS_VERSION_MAJOR
#   TRILINOS_VERSION_MINOR
#   TRILINOS_VERSION_SUBMINOR
#   TRILINOS_WITH_MANDATORY_CXX11
#   TRILINOS_WITH_MPI
#   TRILINOS_SUPPORTS_CPP11
#   TRILINOS_HAS_C99_TR1_WORKAROUND
#

#
# Include the trilinos package configuration:
#
FIND_PACKAGE(Trilinos
  CONFIG
  NAMES Trilinos TRILINOS TrilinosConfig TrilinosConfig.cmake
  HINTS
    ${TRILINOS_DIR}/lib/cmake/Trilinos
    ${TRILINOS_DIR}
    ${CMAKE_SOURCE_DIR}/external/Trilinos/
    $ENV{TRILINOS_DIR}
  PATH_SUFFIXES
    lib64/cmake/Trilinos
    lib/cmake/Trilinos
    lib${LIB_SUFFIX}/cmake/Trilinos
  )

IF(DEFINED Trilinos_VERSION)
  #
  # Extract version numbers:
  #
  SET(TRILINOS_VERSION "${Trilinos_VERSION}")

  STRING(REGEX REPLACE
    "^([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MAJOR "${Trilinos_VERSION}")

  STRING(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MINOR "${Trilinos_VERSION}")

  # If there is no subminor number,
  # TRILINOS_VERSION_SUBMINOR is set to an empty string.
  # If that is the case, set the subminor number to zero
  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.?(([0-9]+)?).*$" "\\1"
    TRILINOS_VERSION_SUBMINOR "${Trilinos_VERSION}")
  IF("${TRILINOS_VERSION_SUBMINOR}" STREQUAL "")
    SET(TRILINOS_VERSION_SUBMINOR "0")
  ENDIF()
ENDIF()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Trilinos DEFAULT_MSG  Trilinos_INCLUDE_DIRS Trilinos_LIBRARIES Trilinos_TPL_LIBRARIES)# Trilinos_TPL_INCLUDE_DIRS Trilinos_additional_headers)

