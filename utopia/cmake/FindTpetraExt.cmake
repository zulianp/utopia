find_path(TPETRAEXT_INCLUDE_DIR
  NAMES TpetraExt_TripleMatrixMultiply.hpp
  HINTS
    ${TRILINOS_DIR}
    ${TPETRAEXT_DIR}
    ${TPETRAEXT_DIR}
    $ENV{TPETRAEXT_DIR}
    $ENV{TRILINOS_DIR}
  PATH_SUFFIXES
    include
  )

if(TPETRAEXT_INCLUDE_DIR)
  set(TPETRAEXT_FOUND TRUE)
else()
  message(WARNING "TPETRAEXT not found")
endif()