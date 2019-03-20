list(APPEND HOMEMADE_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR HOMEMADE_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_utopiafe_module(utopia_fe_homemade utopia_fe_intersection utopia_fe_core)
