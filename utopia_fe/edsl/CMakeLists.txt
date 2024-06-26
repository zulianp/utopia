list(
    APPEND
    EDSL_MODULES
    .
    eval
    eval/mech
    expressions
    local)

# utopiafe_add_module(utopia_edsl ${CMAKE_CURRENT_SOURCE_DIR} "${EDSL_MODULES}")
# target_link_libraries(utopia_edsl PUBLIC Utopia::utopia)

scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${EDSL_MODULES}" UTOPIA_FE_BUILD_INCLUDES UTOPIA_FE_HEADERS UTOPIA_FE_SOURCES)

set(UTOPIA_FE_BUILD_INCLUDES
    ${UTOPIA_FE_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_FE_HEADERS
    ${UTOPIA_FE_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_FE_SOURCES
    ${UTOPIA_FE_SOURCES}
    PARENT_SCOPE)
