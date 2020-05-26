list(APPEND MULTIPHYISCS_MODULES
  .
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR MULTIPHYISCS_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_module(utopia_fe_multiphysics utopia_fe_fractures utopia_fe_core utopia_fe_libmesh_backend transfer utopia_fe_functions utopia_fe_ui utopia_fe_materials)
