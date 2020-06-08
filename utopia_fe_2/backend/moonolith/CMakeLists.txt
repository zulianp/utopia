
if(ENABLE_MOONOLITH)
    list(APPEND INTERSECTION_MODULES
      .
      clipper
      )

    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")
    find_project_files(CMAKE_CURRENT_SOURCE_DIR INTERSECTION_MODULES LOCAL_HEADERS LOCAL_SOURCES)
    #----------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------

    find_package(ParMoonolith)
    #----------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------

    if(MOONOLITH_FOUND)

        add_utopiafe_module(utopia_moonolith utopia_fe_base)
        add_dependencies(utopia_moonolith par_moonolith)
        target_include_directories(utopia_moonolith BEFORE PUBLIC ${MOONOLITH_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/clipper)
        target_link_libraries(utopia_moonolith INTERFACE ${MOONOLITH_LIBRARIES})

        set(MOONOLITH_LIBRARIES ${MOONOLITH_LIBRARIES} PARENT_SCOPE)
        set(MOONOLITH_INCLUDES ${MOONOLITH_INCLUDES} PARENT_SCOPE)
    endif()

endif()
