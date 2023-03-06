list(APPEND UTILS_MODULES . action unit_testing app_management)
# utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${UTILS_MODULES}")

scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${UTILS_MODULES}"
                 UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)

# if(UTOPIA_ENABLE_TRILINOS)
#   # we depend on kokkos headers, i.e. we need to pull in the kokkos target
#   target_link_libraries(utopia PUBLIC Kokkos::kokkos)
# endif()

