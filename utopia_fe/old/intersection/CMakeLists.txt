list(APPEND INTERSECTION_MODULES
  .
  clipper
)
utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${INTERSECTION_MODULES}")

find_package(ParMoonolith REQUIRED)
target_link_libraries(utopia_fe PUBLIC ParMoonolith::par_moonolith)

set(MOONOLITH_LIBRARIES ${MOONOLITH_LIBRARIES} PARENT_SCOPE)
set(MOONOLITH_INCLUDES ${MOONOLITH_INCLUDES} PARENT_SCOPE)
