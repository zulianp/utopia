list(APPEND INTERSECTION_MODULES . clipper)
utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${INTERSECTION_MODULES}")

find_package(ParMoonolith REQUIRED)
target_link_libraries(utopia_fe PUBLIC ParMoonolith::par_moonolith)

set(MOONOLITH_LIBRARIES
    ${MOONOLITH_LIBRARIES}
    PARENT_SCOPE)
set(MOONOLITH_INCLUDES
    ${MOONOLITH_INCLUDES}
    PARENT_SCOPE)

if(ParMoonolith_FOUND)
    list(APPEND MOONOLITH_MODULES .)
    utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${MOONOLITH_MODULES}")
    set(UTOPIA_WITH_MOONOLITH TRUE)
    set(UTOPIA_WITH_MOONOLITH
        TRUE
        PARENT_SCOPE)
endif()
