list(APPEND UTILS_MODULES
    .
    action
    unit_testing
    app_management
)
utopia_add_library(utopia_utils ${CMAKE_CURRENT_SOURCE_DIR} "${UTILS_MODULES}" "")

# almost every module depends on utopia_utils, so we include the base source directory and base build directory
target_include_directories(utopia_utils PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> $<INSTALL_INTERFACE:include>)
target_include_directories(utopia_utils PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>)

if (WITH_TRILINOS)
    # we depend on kokkos headers, i.e. we need to pull in the kokkos target
    target_link_libraries(utopia_utils PUBLIC Kokkos::kokkos)
endif()

find_package(MPIExtended)
if(TRY_WITH_MPI)
    if(MPI_FOUND)
        list(APPEND UTOPIA_INCLUDES ${MPI_CXX_INCLUDE_PATH})
        list(APPEND UTOPIA_LIBRARIES  ${MPI_CXX_LIBRARIES} )

        target_include_directories(utopia_utils SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
        # target_link_libraries(utopia_utils utopia_core ${MPI_CXX_LIBRARIES} )
        target_link_libraries(utopia_utils PUBLIC ${MPI_CXX_LIBRARIES})

        message(STATUS "WHAT: ${MPI_CXX_INCLUDE_PATH} ${MPI_CXX_LIBRARIES}")

        set(WITH_MPI TRUE PARENT_SCOPE)
        if(MPI_DIR)
            set(MPI_DIR ${MPI_DIR} PARENT_SCOPE)
        endif()
    endif()
else()
    message(WARNING "NO Proper MPI installation")
endif()
