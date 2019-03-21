list(APPEND MULTIGRID_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR MULTIGRID_MODULES LOCAL_HEADERS LOCAL_SOURCES)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
add_utopiafe_module(utopia_fe_multigrid transfer utopia_fe_homemade utopia_fe_libmesh_backend utopia_fe_contact)

