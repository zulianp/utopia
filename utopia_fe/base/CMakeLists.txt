list(APPEND BASE_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR BASE_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_module(utopia_fe_base)
