# find dependencies ###################################################

if(BLAS_FOUND OR OPEN_BLAS_FOUND)
  list(APPEND BLAS_MODULES . solvers eigensolvers eval)

  if(LAPACK_FOUND)
    list(APPEND BLAS_MODULES solvers/linear/lapack eigensolvers/lapack)
  endif()

  if(UMFPACK_FOUND)
    list(APPEND BLAS_MODULES solvers/linear/umfpack)
  endif()

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${BLAS_MODULES}"
                   UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

  set(UTOPIA_BUILD_INCLUDES
      ${UTOPIA_BUILD_INCLUDES}
      PARENT_SCOPE)

  set(UTOPIA_HEADERS
      ${UTOPIA_HEADERS}
      PARENT_SCOPE)

  set(UTOPIA_SOURCES
      ${UTOPIA_SOURCES}
      PARENT_SCOPE)

endif()
