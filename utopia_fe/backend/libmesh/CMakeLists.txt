if(LIBMESH_FOUND)
    list(
        APPEND
        UTOPIA_LIBMESH_MODULES
        .
       )

    if(UTOPIA_ENABLE_LIBMESH_LEGACY)
        list(
            APPEND
            UTOPIA_LIBMESH_MODULES
            legacy
           )
    endif()

    if(UTOPIA_ENABLE_LIBMESH_DEPRECATED)
        list(
            APPEND
            UTOPIA_LIBMESH_MODULES
            deprecated
            deprecated/multiphysics
            deprecated/fractures
            deprecated/wear
            deprecated/ui
            deprecated/multigrid
            deprecated/transfer
            deprecated/visualization
            deprecated/intersection
            deprecated/contact
            deprecated/materials
            deprecated/functions)
    endif()



scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${UTOPIA_LIBMESH_MODULES}"
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
endif()

#     utopiafe_add_module(utopia_libmesh ${CMAKE_CURRENT_SOURCE_DIR}
#                         "${UTOPIA_LIBMESH_MODULES}")

#     # to avoid conflict with built-in operator
#     target_compile_definitions(utopia_libmesh PUBLIC PETSC_SKIP_CXX_COMPLEX_FIX)

#     target_include_directories(utopia_libmesh BEFORE
#                                PUBLIC ${LIBMESH_INCLUDE_DIR})

#     target_include_directories(utopia_libmesh BEFORE
#                                PUBLIC ${LIBMESH_INCLUDE_DIRS})

#     target_link_libraries(
#         utopia_libmesh
#         PUBLIC Utopia::utopia utopia_fe_core utopia_edsl utopia_fe_frontend
#                # FIXME remove utopia_interop
#                utopia_interop)

#     target_link_libraries(utopia_libmesh INTERFACE ${LIBMESH_LIBRARIES}) #

#     target_include_directories(
#         utopia_libmesh
#         PRIVATE ${PROJECT_SOURCE_DIR}/backend/libmesh/deprecated/rapidxml)

#     set(LIBMESH_INCLUDE_DIRS
#         ${LIBMESH_INCLUDE_DIRS}
#         PARENT_SCOPE)
#     set(LIBMESH_LIBRARIES
#         ${LIBMESH_LIBRARIES}
#         PARENT_SCOPE)

#     set(UTOPIA_ENABLE_LIBMESH TRUE)
#     set(UTOPIA_ENABLE_LIBMESH
#         TRUE
#         PARENT_SCOPE)

#     set(UTOPIA_FE_LIBRARIES
#         ${UTOPIA_FE_LIBRARIES}
#         PARENT_SCOPE)

# endif()
