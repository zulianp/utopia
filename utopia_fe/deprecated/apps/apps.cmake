if(UTOPIA_ENABLE_LIBMESH)

    set(APPS_MODULES "")
    set(UTOPIA_APPS_DIR ${CMAKE_SOURCE_DIR}/deprecated/apps)

    list(APPEND APPS_MODULES .)

    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")
    find_project_files(${UTOPIA_APPS_DIR} "${APPS_MODULES}" LOCAL_HEADERS
                       LOCAL_SOURCES)
    target_sources(
        utopia_fe_exec
        PRIVATE ${LOCAL_SOURCES}
        PRIVATE ${LOCAL_HEADERS})

    # utopia_link_default_targets(utopia_fe_exec)

    target_include_directories(utopia_fe_exec PRIVATE ${UTOPIA_APPS_DIR})
    target_include_directories(utopia_fe_exec PRIVATE .)
    foreach(MODULE ${APPS_MODULES})
        target_include_directories(utopia_fe_exec
                                   PRIVATE ${UTOPIA_APPS_DIR}/${MODULE})
    endforeach()

endif()
