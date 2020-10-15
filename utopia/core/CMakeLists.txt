list(APPEND CORE_MODULES
    .
    eval
    expressions
    expressions_parallel
    interfaces
    view
    view/expressions
    console
)

utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${CORE_MODULES}")
