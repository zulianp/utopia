list(APPEND TRANSFER_MODULES
  .
  assemblers
  contact_assemblers
  quadrature
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR TRANSFER_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_module(transfer utopia_fe_core utopia_fe_libmesh_backend utopia_fe_intersection utopia_fe_visualization)
foreach(MODULE ${TRANSFER_MODULES})
	target_include_directories(transfer BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)
