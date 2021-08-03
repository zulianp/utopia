list(APPEND CORE_MODULES . fe)

utopiafe_add_module(utopia_fe_core ${CMAKE_CURRENT_SOURCE_DIR}
                    "${CORE_MODULES}")

target_link_libraries(utopia_fe_core PUBLIC Utopia::utopia)
