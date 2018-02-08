list(APPEND MATERIALS_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR MATERIALS_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_materials STATIC ${LOCAL_SOURCES})

target_include_directories(utopia_fe_materials PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(utopia_fe_materials utopia_fe_core utopia_fe_libmesh_backend)


set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_materials" PARENT_SCOPE)

install(TARGETS utopia_fe_materials
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)
