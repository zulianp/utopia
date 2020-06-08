
find_package(LIBMESH)

if(LIBMESH_FOUND)
    list(APPEND LIBMESH_BACKEND_MODULES
        .
        )

    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")
    find_project_files(CMAKE_CURRENT_SOURCE_DIR LIBMESH_BACKEND_MODULES LOCAL_HEADERS LOCAL_SOURCES)

    add_utopiafe_module(utopia_libmesh utopia_fe_base utopia_fe_core)


    target_include_directories(utopia_libmesh BEFORE PUBLIC ${LIBMESH_INCLUDE_DIRS} ${LIBMESH_INCLUDE_DIR})
    target_link_libraries(utopia_libmesh INTERFACE ${LIBMESH_LIBRARIES})

    set(LIBMESH_INCLUDE_DIRS ${LIBMESH_INCLUDE_DIRS} PARENT_SCOPE)
    set(LIBMESH_LIBRARIES ${LIBMESH_LIBRARIES} PARENT_SCOPE)
endif()
