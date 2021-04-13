find_package(Intrepid2)
if(Intrepid2_FOUND)
    utopiafe_add_module(utopia_intrepid2 ${CMAKE_CURRENT_SOURCE_DIR} ".")
    target_link_directories(utopia_intrepid2 PUBLIC ${Intrepid2_LIBRARY_DIRS})

    target_link_libraries(utopia_intrepid2 PUBLIC utopia_fe_core Utopia::utopia)
    target_link_libraries(utopia_intrepid2 PUBLIC ${Intrepid2_LIBRARIES})

    set(UTOPIA_WITH_INTREPID2 TRUE)

    set(UTOPIA_WITH_INTREPID2
        TRUE
        PARENT_SCOPE)
else()
    message(WARNING "could not find Intrepid2")
endif()
