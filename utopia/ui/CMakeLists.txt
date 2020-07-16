set(EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external)

find_path(TINY_EXPR_DIR
    NAMES tinyexpr.h
    HINTS ${EXTERNAL_DIR}/tinyexpr
          ${TINY_EXPR_DIR}
          $ENV{TINY_EXPR_DIR}
          ${INSTALL_DIR}/tinyexpr
          $ENV{INSTALL_DIR}/tinyexpr
)

if(TINY_EXPR_DIR)
    # tinyexpr
    add_library(tinyexpr STATIC ${TINY_EXPR_DIR}/tinyexpr.c)
    target_include_directories(tinyexpr PUBLIC $<BUILD_INTERFACE:${TINY_EXPR_DIR}>)
    set(WITH_TINY_EXPR ON PARENT_SCOPE)
    set(UTOPIA_ADDITIONAL_COMPONENTS ";tinyexpr")

    install(TARGETS tinyexpr
        EXPORT UtopiaTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
endif()

utopia_add_library(utopia_ui ${CMAKE_CURRENT_SOURCE_DIR} "." "utopia_utils;utopia_core")
#hidden dependencies
target_include_directories(utopia_ui SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rapidxml)
# target_compile_options(utopia_ui PRIVATE "-Wno-cast-function-type")

if(WIN32)
    target_include_directories(utopia_ui PRIVATE windows)
endif()

if(TARGET tinyexpr)
    target_link_libraries(utopia_ui PUBLIC tinyexpr)
endif()

find_path(JSON_DIR
    NAMES json.hpp
    HINTS ${EXTERNAL_DIR}/json
          ${JSON_DIR}
          $ENV{JSON_DIR}
          ${INSTALL_DIR}/json
          $ENV{INSTALL_DIR}/json
    PATH_SUFFIXES single_include/nlohmann
)
if(JSON_DIR)
    target_include_directories(utopia_ui SYSTEM PRIVATE ${JSON_DIR})
    set(WITH_JSON ON PARENT_SCOPE)
endif()

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS}${UTOPIA_ADDITIONAL_COMPONENTS}" PARENT_SCOPE)
