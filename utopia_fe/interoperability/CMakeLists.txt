# ##############################################################################
# Interoperability #############################################################
# ##############################################################################

list(APPEND INTEROP_MODULES .)

set(INTEROP_DEP "")
set(INTEROP_EXISTS FALSE)

# ##############################################################################
# ##############################################################################
# Gather dependencies in non-optimal (but not combinatorial) way

if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND INTEROP_DEP utopia_libmesh)
endif()

if(UTOPIA_ENABLE_STK)
    list(APPEND INTEROP_DEP utopia_stk)
endif()

if(UTOPIA_ENABLE_MOONOLITH)
    list(APPEND INTEROP_DEP utopia_moonolith)
endif()

if(UTOPIA_ENABLE_INTREPID2)
    list(APPEND INTEROP_DEP utopia_intrepid2)
endif()

# ##############################################################################
# ##############################################################################

if(UTOPIA_ENABLE_MOONOLITH)

    if(UTOPIA_ENABLE_STK)
        list(APPEND INTEROP_MODULES moonolith_stk)
        set(INTEROP_EXISTS TRUE)
    endif()

    if(UTOPIA_ENABLE_LIBMESH)
        list(
            APPEND
            INTEROP_MODULES
            moonolith_libmesh
            moonolith_libmesh/old
            moonolith_libmesh/old/transfer
            moonolith_libmesh/old/transfer/assemblers
            moonolith_libmesh/old/transfer/quadrature
            moonolith_libmesh/old/transfer/contact_assemblers
            moonolith_libmesh/old/transfer/adapters
            moonolith_libmesh/old/clipper)
        set(INTEROP_EXISTS TRUE)
    endif()

endif()

# ##############################################################################
# ##############################################################################

if(UTOPIA_ENABLE_INTREPID2)
    list(APPEND INTEROP_MODULES intrepid2_generic)

    if(UTOPIA_ENABLE_STK)
        list(APPEND INTEROP_MODULES stk_intrepid2)
        set(INTEROP_EXISTS TRUE)
    endif()

    if(UTOPIA_ENABLE_LIBMESH)
        list(APPEND INTEROP_MODULES libmesh_intrepid2)
    endif()

endif()

# ##############################################################################
# ##############################################################################

if(INTEROP_EXISTS)
    utopiafe_add_module(utopia_interop ${CMAKE_CURRENT_SOURCE_DIR}
                        "${INTEROP_MODULES}")

    message(STATUS "INTEROP_DEP=${INTEROP_DEP}")

    target_link_libraries(utopia_interop PUBLIC Utopia::utopia)
    target_link_libraries(utopia_interop PUBLIC ${INTEROP_DEP})

    set(INTEROP_LIBRARIES
        ${INTEROP_LIBRARIES}
        PARENT_SCOPE)
    set(INTEROP_INCLUDES
        ${INTEROP_INCLUDES}
        PARENT_SCOPE)

endif()

# ##############################################################################
# ##############################################################################
