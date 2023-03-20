# utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")


scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "."
                 UTOPIA_DEP_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

set(UTOPIA_DEP_INCLUDES
    ${UTOPIA_DEP_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)


# ############### SWIG ##################
# cmake_policy(SET "CMP0086" NEW)

