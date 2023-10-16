# set(SLEPC_DIR $ENV{SLEPC_DIR}) message("*** SLEPC_DIR = ${SLEPC_DIR}")

# set(SLEPC_INCLUDE_DIR_A "${SLEPC_DIR}/${PETSC_ARCH}/include")
# set(SLEPC_INCLUDE_DIR "${SLEPC_DIR}/include") set(SLEPC_LIB_DIR
# "${SLEPC_DIR}/${PETSC_ARCH}/lib")

# message("lib dir = ${SLEPC_DIR}/${PETSC_ARCH}/lib")
# find_library(SLEPC_LIB_SLEPC slepc HINTS ${SLEPC_DIR}/lib)

# # message(${SLEPC_LIB_SLEPC})

# if(SLEPC_DIR AND SLEPC_LIB_SLEPC) set(SLEPC_LIBRARIES ${SLEPC_LIB_SLEPC} CACHE
# STRING "SLEPc libraries" FORCE) set(HAVE_SLEPC 1) set(SLEPC_FOUND ON)
# message("-- Found SLEPc: ${SLEPC_LIBRARIES}") set(SLEPC_INCLUDES
# ${SLEPC_INCLUDE_DIR} ${SLEPC_INCLUDE_DIR_A} CACHE STRING "SLEPc include path"
# FORCE) # mark_as_advanced(SLEPC_DIR SLEPC_LIB_SLEPC SLEPC_INCLUDES
# SLEPC_LIBRARIES) else() message(WARNING "SLEPc not found!" ) endif()


if(UTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL)
  set(SLEPC_DIR "${CMAKE_SOURCE_DIR}/../external/petsc")
endif()

if(NOT DEFINED SLEPC_DIR)
  set(SLEPC_DIR $ENV{SLEPC_DIR})
endif()

# Includes ##
# ##############################################################################
if(EXISTS "${SLEPC_DIR}/include" AND EXISTS
                                     "${SLEPC_DIR}/${PETSC_ARCH}/include")
  set(SLEPC_INCLUDES "${SLEPC_DIR}/include"
                     "${SLEPC_DIR}/${PETSC_ARCH}/include")
  set(SLEPC_FOUND TRUE)
else()
  message(SEND_ERROR "SLEPc includes not found")
endif()

# Library ##
# ##############################################################################
if(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
  set(SLEPC_LIBRARIES "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
elseif(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
  set(SLEPC_LIBRARIES "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
elseif(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.dylib")
  set(SLEPC_LIBRARIES "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.dylib")
else()
  message(STATUS "SLEPc library not found")
endif()

# CMake check and done ##
# ##############################################################################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  SLEPc
  "SLEPc could not be found: be sure to set SLEPC_DIR in your environment variables"
  SLEPC_LIBRARIES
  SLEPC_INCLUDES
  SLEPC_DIR)
