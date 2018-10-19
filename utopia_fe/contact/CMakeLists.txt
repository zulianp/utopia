list(APPEND CONTACT_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR CONTACT_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_contact STATIC ${LOCAL_SOURCES})

target_include_directories(utopia_fe_contact BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(utopia_fe_contact utopia_fe_core utopia_fe_libmesh_backend utopia_fe_materials transfer utopia_fe_multigrid utopia_fe_ui)
set_utopia_compiler_features(utopia_fe_contact)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_contact" PARENT_SCOPE)

install(TARGETS utopia_fe_contact
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)
