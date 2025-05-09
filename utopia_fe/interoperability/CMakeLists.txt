# ##############################################################################
# Interoperability #############################################################
# ##############################################################################

list(APPEND INTEROP_MODULES .)

# ##############################################################################
# ##############################################################################
# Gather dependencies in non-optimal (but not combinatorial) way

# if(UTOPIA_ENABLE_LIBMESH) list(APPEND INTEROP_DEP utopia_libmesh) endif()

# if(UTOPIA_ENABLE_STK) list(APPEND INTEROP_DEP utopia_stk) endif()

# if(UTOPIA_ENABLE_MOONOLITH) list(APPEND INTEROP_DEP utopia_moonolith) endif()

# if(UTOPIA_ENABLE_INTREPID2) list(APPEND INTEROP_DEP utopia_intrepid2) endif()

# if(UTOPIA_ENABLE_PETSC) list(APPEND INTEROP_DEP utopia_petsc_dm) endif()

# ##############################################################################
# ##############################################################################

if(UTOPIA_ENABLE_MOONOLITH)

  list(APPEND INTEROP_MODULES moonolith_generic)

  if(UTOPIA_ENABLE_STK)
    list(APPEND INTEROP_MODULES moonolith_stk)
  endif()

  if(UTOPIA_ENABLE_PETSCDM)
    list(APPEND INTEROP_MODULES moonolith_petsc)
  endif()

  if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND INTEROP_MODULES moonolith_libmesh)

    if(UTOPIA_ENABLE_LIBMESH_LEGACY)
      list(APPEND INTEROP_MODULES moonolith_libmesh/legacy)
    endif()

    if(UTOPIA_ENABLE_LIBMESH_DEPRECATED)
      list(
        APPEND
        INTEROP_MODULES
        moonolith_libmesh/deprecated
        moonolith_libmesh/deprecated/transfer
        moonolith_libmesh/deprecated/transfer/assemblers
        moonolith_libmesh/deprecated/transfer/quadrature
        moonolith_libmesh/deprecated/transfer/contact_assemblers
        moonolith_libmesh/deprecated/transfer/adapters
        moonolith_libmesh/deprecated/clipper
        moonolith_libmesh/deprecated/contact
        moonolith_libmesh/deprecated/transfer
        # moonolith_libmesh/deprecated/multigrid
        # moonolith_libmesh/deprecated/wear
      )
    endif()

  endif()

endif()

# ##############################################################################
# ##############################################################################

if(UTOPIA_ENABLE_INTREPID2)
  list(APPEND INTEROP_MODULES intrepid2_generic)

  if(UTOPIA_ENABLE_STK)
    list(APPEND INTEROP_MODULES stk_intrepid2)

  endif()

  if(UTOPIA_ENABLE_ARBORX)
    list(APPEND INTEROP_MODULES intrepid2_arborx)

  endif()

  if(UTOPIA_ENABLE_PETSCDM)
    list(APPEND INTEROP_MODULES petsc_intrepid2)

  endif()

endif()

if(UTOPIA_ENABLE_SFEM)
  list(APPEND INTEROP_MODULES sfem_generic)

  if(UTOPIA_ENABLE_STK AND UTOPIA_ENABLE_MOONOLITH)
    list(APPEND INTEROP_MODULES sfem_stk)

    list(APPEND INTEROP_MODULES utopia_moonolith)
  endif()
endif()

if(UTOPIA_ENABLE_LIBMESH
   AND UTOPIA_ENABLE_KOKKOS
   AND UTOPIA_ENABLE_LIBMESH_KOKKOS)
  list(APPEND INTEROP_MODULES libmesh_kokkos)
endif()

# ##############################################################################
# ##############################################################################

# utopiafe_add_module(utopia_interop ${CMAKE_CURRENT_SOURCE_DIR}
# "${INTEROP_MODULES}")

# message(STATUS "INTEROP_DEP=${INTEROP_DEP}")

# target_link_libraries(utopia_interop PUBLIC Utopia::utopia)
# target_link_libraries(utopia_interop PUBLIC ${INTEROP_DEP})

# set(INTEROP_LIBRARIES ${INTEROP_LIBRARIES} PARENT_SCOPE)

# set(INTEROP_INCLUDES ${INTEROP_INCLUDES} PARENT_SCOPE)

# set(UTOPIA_FE_LIBRARIES ${UTOPIA_FE_LIBRARIES} PARENT_SCOPE)

scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${INTEROP_MODULES}"
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

# ##############################################################################
# ##############################################################################
