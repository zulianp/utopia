utopiafe_add_module(utopia_petsc_dm ${CMAKE_CURRENT_SOURCE_DIR} ".")
target_link_libraries(utopia_petsc_dm PUBLIC utopia_fe_core utopia_fe_frontend Utopia::utopia)

set(UTOPIA_WITH_PETSCDM TRUE)

set(UTOPIA_WITH_PETSCDM
    TRUE
    PARENT_SCOPE)

set(UTOPIA_FE_LIBRARIES
    ${UTOPIA_FE_LIBRARIES}
    PARENT_SCOPE)