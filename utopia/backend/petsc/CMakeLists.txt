# set-up module
if(PETSC_FOUND)
  list(
    APPEND
    PETSC_MODULES
    .
    eval
    operators
    solvers
    solvers/linear
    solvers/multilevel
    solvers/nonlinear
    solvers/smoothers
    solvers/nonlinear/constrained
    solvers/nonlinear/constrained/quadratic_programming
    solvers/nonlinear/trust_region
    solvers/dd
    views)

  if(UTOPIA_ENABLE_SLEPC)
    list(APPEND PETSC_MODULES solvers/eigensolvers)
  endif()

  if(UTOPIA_ENABLE_PASSO_EXTENSIONS)
    list(APPEND PETSC_MODULES external)
  endif()

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${PETSC_MODULES}"
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
