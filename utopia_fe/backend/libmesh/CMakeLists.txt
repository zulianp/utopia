
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
