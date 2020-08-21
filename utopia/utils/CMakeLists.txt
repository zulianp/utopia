list(APPEND UTILS_MODULES
    .
    action
    unit_testing
    app_management
)
utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${UTILS_MODULES}")

if (UTOPIA_WITH_TRILINOS)
    # we depend on kokkos headers, i.e. we need to pull in the kokkos target
    target_link_libraries(utopia PUBLIC Kokkos::kokkos)
endif()

find_package(MPIExtended)
if(UTOPIA_ENABLE_MPI)
    if(MPI_FOUND)
        list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${MPI_CXX_INCLUDE_PATH})
        list(APPEND UTOPIA_THIRDPARTY_LIBRARIES  ${MPI_CXX_LIBRARIES} )

        target_include_directories(utopia SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
        # target_link_libraries(utopia_utils utopia_core ${MPI_CXX_LIBRARIES} )
        target_link_libraries(utopia PUBLIC ${MPI_CXX_LIBRARIES})

        message(STATUS "WHAT: ${MPI_CXX_INCLUDE_PATH} ${MPI_CXX_LIBRARIES}")

        set(UTOPIA_WITH_MPI TRUE PARENT_SCOPE)
        if(MPI_DIR)
            set(MPI_DIR ${MPI_DIR} PARENT_SCOPE)
        endif()
    endif()
else()
    message(WARNING "NO Proper MPI installation")
endif()
