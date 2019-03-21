list(APPEND CORE_MODULES
	.
	eval
	expressions
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR CORE_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_module(utopia_fe_core)
foreach(MODULE ${CORE_MODULES})
	target_include_directories(utopia_fe_core BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)
