list(APPEND UI_MODULES
    .
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR UI_MODULES LOCAL_HEADERS LOCAL_SOURCES)

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
    target_include_directories(tinyexpr PUBLIC ${TINY_EXPR_DIR})
    set(WITH_TINY_EXPR ON PARENT_SCOPE)
    set(UTOPIA_ADDITIONAL_COMPONENTS ";tinyexpr")

    install(TARGETS tinyexpr
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
endif()


add_library(utopia_ui STATIC ${LOCAL_SOURCES})
utopia_link_default_targets(utopia_ui)

target_include_directories(utopia_ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})
target_include_directories(utopia_ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(utopia_ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

#hidden dependencies
target_include_directories(utopia_ui SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rapidxml)
# target_include_directories(utopia_ui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rapidjson)

if(WIN32)
    target_include_directories(utopia_ui PRIVATE windows)
endif()

target_link_libraries(utopia_ui utopia_utils utopia_core)

if(TARGET tinyexpr)
    target_link_libraries(utopia_ui tinyexpr)
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

set_utopia_compiler_features(utopia_ui)


install(TARGETS utopia_ui
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_ui${UTOPIA_ADDITIONAL_COMPONENTS}" PARENT_SCOPE)

# message(STATUS "${UTOPIA_COMPONENTS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")


set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
