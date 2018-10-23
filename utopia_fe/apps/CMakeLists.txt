list(APPEND APPS_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR APPS_MODULES LOCAL_HEADERS LOCAL_SOURCES)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_core" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")


add_library(utopia_fe_apps STATIC ${LOCAL_SOURCES})
target_link_libraries(
	utopia_fe_apps
	transfer
	utopia_fe_visualization
	utopia_fe_homemade
	utopia_fe_libmesh_backend
	utopia_fe_contact
	utopia_fe_materials
	utopia_fe_multigrid
	utopia_fe_functions
	utopia_fe_fractures
	utopia_fe_wear
)

target_include_directories(utopia_fe_apps BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/simulations)
set_utopia_compiler_features(utopia_fe_apps)

install(TARGETS utopia_fe_apps
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_apps" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
