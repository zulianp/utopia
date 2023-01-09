list(
    APPEND
    CORE_MODULES
    .
    functions
    eval
    expressions
    expressions_parallel
    interfaces
    view
    view/expressions
    console
    misc
    traits)

# utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${CORE_MODULES}")
scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${CORE_MODULES}" UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)


set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)
