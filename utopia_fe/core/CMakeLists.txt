list(APPEND CORE_MODULES
    .
    eval
    eval/mech
    expressions
    local
)
utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${CORE_MODULES}")
