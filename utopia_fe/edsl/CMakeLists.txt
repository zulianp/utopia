list(
    APPEND
    EDSL_MODULES
    .
    eval
    eval/mech
    expressions
    local)

utopiafe_add_module(utopia_edsl ${CMAKE_CURRENT_SOURCE_DIR} "${EDSL_MODULES}")
target_link_libraries(utopia_edsl PUBLIC Utopia::utopia)
