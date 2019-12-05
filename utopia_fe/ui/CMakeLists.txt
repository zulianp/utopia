list(APPEND UI_MODULES
  .
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR UI_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_module(utopia_fe_ui utopia_fe_core utopia_fe_materials utopia_fe_contact par_moonolith)


foreach(MODULE ${UI_MODULES})
    target_include_directories(utopia_fe_ui BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)
