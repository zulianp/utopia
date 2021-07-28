if(UTOPIA_ENABLE_ARBORX)
    find_package(ArborX REQUIRED)
    utopiafe_add_module(utopia_arborx ${CMAKE_CURRENT_SOURCE_DIR} ".")
    target_link_libraries(utopia_arborx ArborX::ArborX)


endif()
