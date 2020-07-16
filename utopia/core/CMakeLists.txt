list(APPEND CORE_MODULES
    .
    eval
    expressions
    expressions_parallel
    interfaces
    view
    view/expressions
)

utopia_add_library(utopia_core ${CMAKE_CURRENT_SOURCE_DIR} "${CORE_MODULES}" "utopia_utils;utopia_ui" "additional_file.cpp;add2.cpp")
