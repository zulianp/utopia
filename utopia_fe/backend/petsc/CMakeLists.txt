# if(PETSc_FOUND)
  list(APPEND PETSC_MODULES .)

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${PETSC_MODULES}"
                   UTOPIA_FE_BUILD_INCLUDES UTOPIA_FE_HEADERS UTOPIA_FE_SOURCES)

  set(UTOPIA_FE_BUILD_INCLUDES
      ${UTOPIA_FE_BUILD_INCLUDES}
      PARENT_SCOPE)

  set(UTOPIA_FE_HEADERS
      ${UTOPIA_FE_HEADERS}
      PARENT_SCOPE)

  set(UTOPIA_FE_SOURCES
      ${UTOPIA_FE_SOURCES}
      PARENT_SCOPE)
# endif()


# utopiafe_add_module(utopia_petsc_dm ${CMAKE_CURRENT_SOURCE_DIR} ".")
# target_link_libraries(utopia_petsc_dm PUBLIC utopia_fe_core utopia_fe_frontend
#                                              Utopia::utopia)

# set(UTOPIA_ENABLE_PETSCDM TRUE)

# set(UTOPIA_ENABLE_PETSCDM
#     TRUE
#     PARENT_SCOPE)

# set(UTOPIA_FE_LIBRARIES
#     ${UTOPIA_FE_LIBRARIES}
#     PARENT_SCOPE)
