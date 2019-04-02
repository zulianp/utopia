list(APPEND TEST_MODULES
    .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR TEST_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_module(utopia_fe_tests
    transfer
    utopia_fe_visualization
    utopia_fe_homemade
    utopia_fe_intrepid2_backend
    utopia_fe_libmesh_backend
    utopia_fe_contact
    utopia_fe_materials
    utopia_fe_multigrid
    utopia_fe_functions
    utopia_fe_fractures
    utopia_fe_wear
    utopia_fe_apps
)

if(TARGET utopia_fe_intrepid2)
    message(STATUS "compiling intrepid2 tests")
    target_link_libraries(utopia_fe_tests INTERFACE utopia_fe_intrepid2_backend)
endif()
