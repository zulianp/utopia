# find_path(TRILINOS_TPETRAEXT_INCLUDE_DIR
#   NAMES TpetraExt_TripleMatrixMultiply.hpp
#   HINTS
#     ${TRILINOS_DIR}
#     ${TRILINOS_TPETRAEXT_DIR}
#     $ENV{TRILINOS_TPETRAEXT_DIR}
#     $ENV{TRILINOS_DIR}
#   PATH_SUFFIXES
#     include
#   )

# if(TRILINOS_TPETRAEXT_INCLUDE_DIR)
#   set(TRILINOS_TPETRAEXT_FOUND TRUE)
# else()
#   message(WARNING "TRILINOS_TPETRAEXT not found")
# endif()
