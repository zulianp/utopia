list(APPEND WEAR_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR WEAR_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_module(utopia_fe_wear
    utopia_fe_visualization
    # utopia_fe_homemade
    utopia_fe_libmesh_backend
    utopia_fe_contact
    utopia_fe_materials
    utopia_fe_multigrid
    utopia_fe_functions
    utopia_fe_fractures
    utopia_fe_ui
)
