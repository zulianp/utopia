list(APPEND TRANSFER_MODULES
  .
  adapters
  assemblers
  contact_assemblers
  quadrature
)
utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${TRANSFER_MODULES}")
