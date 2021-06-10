list(APPEND FRONTEND_MODULES .)

utopiafe_add_module(utopia_fe_frontend ${CMAKE_CURRENT_SOURCE_DIR}
                    "${FRONTEND_MODULES}")

target_link_libraries(utopia_fe_frontend PUBLIC Utopia::utopia)
