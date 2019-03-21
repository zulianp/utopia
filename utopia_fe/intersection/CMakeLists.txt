list(APPEND INTERSECTION_MODULES
  .
  clipper
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR INTERSECTION_MODULES LOCAL_HEADERS LOCAL_SOURCES)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(ParMoonolith REQUIRED)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

add_utopiafe_module(utopia_fe_intersection utopia_fe_visualization)
add_dependencies(utopia_fe_intersection par_moonolith)
target_include_directories(utopia_fe_intersection BEFORE PUBLIC ${MOONOLITH_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/clipper)
target_link_libraries(utopia_fe_intersection INTERFACE ${MOONOLITH_LIBRARIES})

set(MOONOLITH_LIBRARIES ${MOONOLITH_LIBRARIES} PARENT_SCOPE)
set(MOONOLITH_INCLUDES ${MOONOLITH_INCLUDES} PARENT_SCOPE)
