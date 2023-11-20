list(APPEND FRONTEND_MODULES .
    fracnetflow
    functions
    integrators
    integrators/variants
    multigrid
    solve
    fsi
    mech
    domain_decomposition
    # franetg
)

# utopiafe_add_module(utopia_fe_frontend ${CMAKE_CURRENT_SOURCE_DIR}
#                     "${FRONTEND_MODULES}")

# target_link_libraries(utopia_fe_frontend PUBLIC Utopia::utopia utopia_fe_core)


scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${FRONTEND_MODULES}" UTOPIA_FE_BUILD_INCLUDES UTOPIA_FE_HEADERS UTOPIA_FE_SOURCES)


set(UTOPIA_FE_BUILD_INCLUDES
    ${UTOPIA_FE_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_FE_HEADERS
    ${UTOPIA_FE_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_FE_SOURCES
    ${UTOPIA_FE_SOURCES}
    PARENT_SCOPE)