if(UTOPIA_ENABLE_STK)
    find_package(Intrepid2)
    if(Intrepid2_FOUND)
        utopiafe_add_module(utopia_intrepid2 ${CMAKE_CURRENT_SOURCE_DIR} ".")

        target_include_directories(utopia_intrepid2 BEFORE
                                   PUBLIC ${CMAKE_SOURCE_DIR}/base)

        target_include_directories(utopia_intrepid2 BEFORE
                                   PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

        target_include_directories(utopia_intrepid2 BEFORE
                                   PUBLIC ${CMAKE_BINARY_DIR})

        target_link_libraries(utopia_intrepid2 PUBLIC intrepid2 Utopia::utopia)

        set(UTOPIA_WITH_INTREPID2 TRUE)

        set(UTOPIA_WITH_INTREPID2
            TRUE
            PARENT_SCOPE)
    else()
        message(WARNING "could not find Intrepid2")
    endif()
endif()
