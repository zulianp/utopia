# set-up module
if(Trilinos_FOUND)
  list(
    APPEND
    TRILINOS_MODULES
    .
    eval
    solvers
    solvers/linear
    solvers/nonlinear
    kokkos)

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${TRILINOS_MODULES}"
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

  set(UTOPIA_DEP_INCLUDES
      ${Trilinos_INCLUDE_DIRS}
      PARENT_SCOPE)
endif()
