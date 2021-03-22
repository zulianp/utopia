find_package(ParMoonolith REQUIRED)

if(ParMoonolith_FOUND)
    list(APPEND MOONOLITH_MODULES .)
    utopiafe_add_module(utopia_moonolith ${CMAKE_CURRENT_SOURCE_DIR}
                        "${MOONOLITH_MODULES}")
    set(UTOPIA_WITH_MOONOLITH TRUE)
    set(UTOPIA_WITH_MOONOLITH
        TRUE
        PARENT_SCOPE)

    install(
        TARGETS utopia_moonolith
        EXPORT UtopiaMoonolithTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

    target_link_libraries(utopia_moonolith PUBLIC ParMoonolith::par_moonolith
                                                  Utopia::utopia)

    target_include_directories(utopia_moonolith BEFORE
                               PUBLIC ${CMAKE_SOURCE_DIR}/base)

    target_include_directories(utopia_moonolith BEFORE
                               PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    target_include_directories(utopia_moonolith BEFORE
                               PUBLIC ${CMAKE_BINARY_DIR})

    set(MOONOLITH_LIBRARIES
        ${MOONOLITH_LIBRARIES}
        PARENT_SCOPE)
    set(MOONOLITH_INCLUDES
        ${MOONOLITH_INCLUDES}
        PARENT_SCOPE)
endif()
