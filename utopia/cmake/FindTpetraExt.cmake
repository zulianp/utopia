set(TPETRA_SEARCH_PATHS
    "${Trilinos_DIR}/../../../;${Trilinos_TPETRAEXT_DIR};${Trilinos_DIR};${TRILINOS_INSTALL_DIR}"
)

if(UTOPIA_ENABLE_ENV_READ)
  set(TPETRA_SEARCH_PATHS
      "${TPETRA_SEARCH_PATHS};$ENV{Trilinos_TPETRAEXT_DIR};$ENV{Trilinos_DIR}")

endif()

find_path(
  TRILINOS_TPETRAEXT_INCLUDE_DIR
  NAMES TpetraExt_TripleMatrixMultiply.hpp
  PATHS ${TPETRA_SEARCH_PATHS}
  PATH_SUFFIXES include)

if(TRILINOS_TPETRAEXT_INCLUDE_DIR)
  set(TRILINOS_TPETRAEXT_FOUND TRUE)
else()
  # message(STATUS "Trilinos dir: ${TRILINOS_DIR}")
  message(WARNING "TRILINOS_TPETRAEXT not found")
endif()
