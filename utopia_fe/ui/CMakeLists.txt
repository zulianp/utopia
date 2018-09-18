list(APPEND UI_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR UI_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_ui STATIC ${LOCAL_SOURCES})

target_include_directories(utopia_fe_ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
	utopia_fe_ui utopia_fe_core utopia_fe_materials
)

set_utopia_compiler_features(utopia_fe_ui)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_ui" PARENT_SCOPE)

install(TARGETS utopia_fe_ui
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)
