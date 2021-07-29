if(UTOPIA_ENABLE_ARBORX)
    # find_package(MPI COMPONENTS CXX REQUIRED)

    # FIXME
    set(MPI_C_FOUND TRUE)
    set(MPI_C_WORKS TRUE)
    set(MPI_CXX_WORKS TRUE)
    set(MPI_CXX_FOUND TRUE)

    find_package(ArborX REQUIRED)
    utopiafe_add_module(utopia_arborx ${CMAKE_CURRENT_SOURCE_DIR} ".")
    target_link_libraries(utopia_arborx ArborX::ArborX Utopia::utopia)

endif()
