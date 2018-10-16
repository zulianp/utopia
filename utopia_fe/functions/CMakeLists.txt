list(APPEND FUNCTIONS_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR FUNCTIONS_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_functions STATIC ${LOCAL_SOURCES})

target_include_directories(utopia_fe_functions BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(utopia_fe_functions utopia_fe_core utopia_fe_libmesh_backend)
set_utopia_compiler_features(utopia_fe_functions)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_functions" PARENT_SCOPE)

install(TARGETS utopia_fe_functions
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)
