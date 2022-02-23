find_package(LIBMESH)

if(LIBMESH_FOUND)

    list(
        APPEND
        UTOPIA_LIBMESH_MODULES
        .
        old
        old/multiphysics
        old/fractures
        old/wear
        old/ui
        old/multigrid
        old/transfer
        old/visualization
        backend/libmesh
        old/intersection
        old/contact
        old/materials
        old/functions)

    utopiafe_add_module(utopia_libmesh ${CMAKE_CURRENT_SOURCE_DIR}
                        "${UTOPIA_LIBMESH_MODULES}")

    target_include_directories(utopia_libmesh BEFORE
                               PUBLIC ${LIBMESH_INCLUDE_DIR})

    target_include_directories(utopia_libmesh BEFORE
                               PUBLIC ${LIBMESH_INCLUDE_DIRS})

    target_link_libraries(
        utopia_libmesh
        PUBLIC Utopia::utopia utopia_fe_core utopia_edsl utopia_fe_frontend
               # FIXME remove utopia_interop
               utopia_interop)

    target_link_libraries(utopia_libmesh INTERFACE ${LIBMESH_LIBRARIES}) #

    target_include_directories(
        utopia_libmesh
        PRIVATE ${PROJECT_SOURCE_DIR}/backend/libmesh/old/rapidxml)

    set(LIBMESH_INCLUDE_DIRS
        ${LIBMESH_INCLUDE_DIRS}
        PARENT_SCOPE)
    set(LIBMESH_LIBRARIES
        ${LIBMESH_LIBRARIES}
        PARENT_SCOPE)

    set(UTOPIA_WITH_LIBMESH TRUE)
    set(UTOPIA_WITH_LIBMESH
        TRUE
        PARENT_SCOPE)

    set(UTOPIA_FE_LIBRARIES
        ${UTOPIA_FE_LIBRARIES}
        PARENT_SCOPE)

endif()
