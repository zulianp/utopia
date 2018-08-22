list(APPEND WEAR_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR WEAR_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_wear STATIC ${LOCAL_SOURCES})

target_include_directories(utopia_fe_wear PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
	utopia_fe_wear
	utopia_fe_visualization
	utopia_fe_homemade
	utopia_fe_libmesh_backend
	utopia_fe_contact
	utopia_fe_materials
	utopia_fe_multigrid
	utopia_fe_functions
	utopia_fe_fractures
)

set_utopia_compiler_features(utopia_fe_wear)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_wear" PARENT_SCOPE)

install(TARGETS utopia_fe_wear
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)
