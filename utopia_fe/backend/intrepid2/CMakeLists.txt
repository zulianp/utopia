find_package(Intrepid2)
if(INTREPID2_FOUND)

    list(APPEND INTREPID2_BACKEND_MODULES
        .
    )

    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")
    find_project_files(CMAKE_CURRENT_SOURCE_DIR INTREPID2_BACKEND_MODULES LOCAL_HEADERS LOCAL_SOURCES)

    add_utopiafe_module(utopia_fe_intrepid2_backend utopia_fe_core utopia_fe_intersection utopia_fe_base)
    target_include_directories(utopia_fe_intrepid2_backend PUBLIC ${INTREPID2_INCLUDE_DIR})
    target_link_libraries(utopia_fe_intrepid2_backend INTERFACE ${INTREPID2_LIBRARIES})

    set(WITH_INTREPID2 TRUE PARENT_SCOPE)
else()
    message(WARNING "could not find Intrepid2")
endif()

