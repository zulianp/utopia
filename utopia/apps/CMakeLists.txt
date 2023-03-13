# set(APPS_MODULES "")
# set(UTOPIA_APPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/apps)

# if(UTOPIA_ENABLE_PETSC)
#   list(
#     APPEND
#     APPS_MODULES
#     .
#     data_structures
#     gpu_assembler
#     cpu_assembler
#     cpu_assembler/dmplex
#     PF_frac
#     PF_frac/vc
#     min_problems
#     fe
#     ls_solve
#     qp_solve)

#   if(UTOPIA_ENABLE_VC)
#     # Requires petsc
#     list(APPEND APPS_MODULES simd_assembler)
#     list(APPEND APPS_MODULES simd_assembler_v2)
#   endif()

# endif()


# scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${APPS_MODULES}"
#                  UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

# set(UTOPIA_BUILD_INCLUDES
#     ${UTOPIA_BUILD_INCLUDES}
#     PARENT_SCOPE)

# set(UTOPIA_HEADERS
#     ${UTOPIA_HEADERS}
#     PARENT_SCOPE)

# set(UTOPIA_SOURCES
#     ${UTOPIA_SOURCES}
#     PARENT_SCOPE)