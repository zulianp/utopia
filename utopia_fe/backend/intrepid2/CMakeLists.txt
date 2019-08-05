
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(Intrepid2)
if(INTREPID2_FOUND)

    list(APPEND INTREPID2_BACKEND_MODULES
        .
    )

    find_project_files(CMAKE_CURRENT_SOURCE_DIR INTREPID2_BACKEND_MODULES LOCAL_HEADERS LOCAL_SOURCES)

    add_utopiafe_module(utopia_fe_intrepid2 utopia_fe_core utopia_fe_intersection utopia_fe_base)
    target_include_directories(utopia_fe_intrepid2 PUBLIC ${INTREPID2_INCLUDE_DIR})
    target_link_libraries(utopia_fe_intrepid2 INTERFACE ${INTREPID2_LIBRARIES})

    # set(INTREPID2_INCLUDE_DIR ${Intrepid2_INCLUDE_DIR} PARENT_SCOPE)
    # set(Intrepid2_LIBRARIES ${Intrepid2_LIBRARIES} PARENT_SCOPE)

    #----------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------

    set(WITH_INTREPID2 TRUE PARENT_SCOPE)
else()
    message(WARNING "could not find Intrepid2")
endif()



