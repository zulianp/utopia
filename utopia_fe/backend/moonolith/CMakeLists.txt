find_package(ParMoonolith REQUIRED)

# if(ParMoonolith_FOUND)
list(APPEND MOONOLITH_MODULES .)
utopiafe_add_module(utopia_moonolith ${CMAKE_CURRENT_SOURCE_DIR}
                    "${MOONOLITH_MODULES}")
set(UTOPIA_ENABLE_MOONOLITH TRUE)
set(UTOPIA_ENABLE_MOONOLITH
    TRUE
    PARENT_SCOPE)

target_link_libraries(utopia_moonolith PUBLIC ParMoonolith::par_moonolith
                                              Utopia::utopia utopia_fe_core)

set(MOONOLITH_LIBRARIES
    ${MOONOLITH_LIBRARIES}
    PARENT_SCOPE)
set(MOONOLITH_INCLUDES
    ${MOONOLITH_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_FE_LIBRARIES
    ${UTOPIA_FE_LIBRARIES}
    PARENT_SCOPE)
