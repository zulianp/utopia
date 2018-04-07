set(SLEPC_DIR $ENV{SLEPC_DIR})
message("*** SLEPC_DIR = ${SLEPC_DIR}")


set(SLEPC_INCLUDE_DIR_A "${SLEPC_DIR}/${PETSC_ARCH}/include")
set(SLEPC_INCLUDE_DIR   "${SLEPC_DIR}/include")
set(SLEPC_LIB_DIR       "${SLEPC_DIR}/${PETSC_ARCH}/lib")

message(" lib dir = ${SLEPC_DIR}/${PETSC_ARCH}/lib")
find_library(SLEPC_LIB_SLEPC     slepc 
             HINTS ${SLEPC_DIR}/${PETSC_ARCH}/lib)

message(${SLEPC_LIB_SLEPC})

if (SLEPC_DIR AND SLEPC_LIB_SLEPC )
  set(SLEPC_LIBRARIES ${SLEPC_LIB_SLEPC} CACHE STRING "SLEPc libraries" FORCE)
  set(HAVE_SLEPC 1)
  set(SLEPC_FOUND ON)
  message( "-- Found SLEPc: ${SLEPC_LIBRARIES}" )
  set(SLEPC_INCLUDES ${SLEPC_INCLUDE_DIR} ${SLEPC_INCLUDE_DIR_A} CACHE STRING "SLEPc include path" FORCE)
  mark_as_advanced( SLEPC_DIR SLEPC_LIB_SLEPC SLEPC_INCLUDES SLEPC_LIBRARIES )
else()
  message(WARNING "SLEPc not found!" )
endif()
