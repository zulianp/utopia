list(
    APPEND
    CORE_MODULES
    .
    functions
    eval
    expressions
    expressions_parallel
    interfaces
    view
    view/expressions
    console
    misc
    traits)

utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${CORE_MODULES}")
