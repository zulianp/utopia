set(EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external)

# ##############################################################################
# TINYEXPR

if(UTOPIA_ENABLE_TINY_EXPR)

  find_path(
    TINY_EXPR_DIR
    NAMES tinyexpr.h
    HINTS ${EXTERNAL_DIR}/tinyexpr ${TINY_EXPR_DIR} $ENV{TINY_EXPR_DIR}
          ${INSTALL_DIR}/tinyexpr $ENV{INSTALL_DIR}/tinyexpr)

  if(TINY_EXPR_DIR)
    # tinyexpr
    add_library(tinyexpr ${TINY_EXPR_DIR}/tinyexpr.c)
    target_include_directories(tinyexpr
                               PUBLIC $<BUILD_INTERFACE:${TINY_EXPR_DIR}>)


    # Add headers and sources to global variables.
    scan_directories(${TINY_EXPR_DIR} "." UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS
                     UTOPIA_SOURCES)
    set(UTOPIA_BUILD_INCLUDES
        ${UTOPIA_BUILD_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_HEADERS
        ${UTOPIA_HEADERS}
        PARENT_SCOPE)

    set(UTOPIA_SOURCES
        ${UTOPIA_SOURCES}
        PARENT_SCOPE)

    set(UTOPIA_ENABLE_TINY_EXPR
        ON
        PARENT_SCOPE)
    set(UTOPIA_ADDITIONAL_COMPONENTS ";tinyexpr")

    install(
      TARGETS tinyexpr
      EXPORT UtopiaTargets
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
  endif()
endif()




# TODO: ASK if rapidxml should be optional.


# ##############################################################################

# utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")


list(APPEND 
    UI_MOUDULES
    .
    rapidjson
    rapidxml
    )

scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${UI_MOUDULES}" UTOPIA_BUILD_INCLUDES
                 UTOPIA_HEADERS UTOPIA_SOURCES)


# hidden dependencies target_include_directories(utopia SYSTEM PRIVATE
# ${CMAKE_CURRENT_SOURCE_DIR}/rapidxml)

# ##############################################################################
# target_compile_options(utopia PRIVATE "-Wno-cast-function-type")

if(WIN32)
  # target_include_directories(utopia PRIVATE windows)
  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "windows" UTOPIA_BUILD_INCLUDES
                 UTOPIA_HEADERS UTOPIA_SOURCES)

endif()

# TODO REMEMBER TO FIX if(TARGET tinyexpr) target_link_libraries(utopia PUBLIC
# tinyexpr) endif()

# SAME if(UTOPIA_ENABLE_RAPIDJSON) find_path( JSON_DIR NAMES json.hpp HINTS
# ${EXTERNAL_DIR}/json ${JSON_DIR} $ENV{JSON_DIR} ${INSTALL_DIR}/json
# $ENV{INSTALL_DIR}/json PATH_SUFFIXES single_include/nlohmann) if(JSON_DIR)
# target_include_directories(utopia SYSTEM PRIVATE ${JSON_DIR})
# set(UTOPIA_ENABLE_JSON ON PARENT_SCOPE) endif() endif()

set(UTOPIA_COMPONENTS
    "${UTOPIA_COMPONENTS}${UTOPIA_ADDITIONAL_COMPONENTS}"
    PARENT_SCOPE)


set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)
