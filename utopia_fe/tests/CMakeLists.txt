list(APPEND TEST_MODULES
    .
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(${CMAKE_CURRENT_SOURCE_DIR} "${TEST_MODULES}" LOCAL_HEADERS LOCAL_SOURCES)
add_utopiafe_app(utopia_fe_tests
    utopia_fe
    utopia_fe_apps
)
