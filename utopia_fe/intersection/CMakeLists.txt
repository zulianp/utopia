list(APPEND INTERSECTION_MODULES
  .
  clipper
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR INTERSECTION_MODULES LOCAL_HEADERS LOCAL_SOURCES)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(ParMoonolith)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

add_library(utopia_fe_intersection STATIC ${LOCAL_SOURCES})
add_dependencies(utopia_fe_intersection par_moonolith)
target_include_directories(utopia_fe_intersection BEFORE PUBLIC ${MOONOLITH_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/clipper)
target_link_libraries(utopia_fe_intersection ${MOONOLITH_LIBRARIES})
target_link_libraries(utopia_fe_intersection ${MOONOLITH_LIBRARIES})
target_link_libraries(utopia_fe_intersection utopia_fe_visualization)
set_utopia_compiler_features(utopia_fe_intersection)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_intersection" PARENT_SCOPE)

set(MOONOLITH_LIBRARIES ${MOONOLITH_LIBRARIES} PARENT_SCOPE)
set(MOONOLITH_INCLUDES ${MOONOLITH_INCLUDES} PARENT_SCOPE)

install(TARGETS utopia_fe_intersection
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)
