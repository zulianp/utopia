list(APPEND FRACTURES_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR FRACTURES_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_fractures STATIC ${LOCAL_SOURCES})

target_include_directories(utopia_fe_fractures BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(utopia_fe_fractures utopia_fe_core utopia_fe_libmesh_backend transfer)
set_utopia_compiler_features(utopia_fe_fractures)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_fractures" PARENT_SCOPE)

install(TARGETS utopia_fe_fractures
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)
