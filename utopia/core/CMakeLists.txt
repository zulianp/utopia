list(APPEND CORE_MODULES
    .
    eval
    expressions
    expressions_parallel
    interfaces
    view
    view/expressions
    console
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR CORE_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_core STATIC utopia_core.cpp ${LOCAL_SOURCES})
utopia_link_default_targets(utopia_core)

target_include_directories(utopia_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})
foreach(MODULE ${CORE_MODULES})
    target_include_directories(utopia_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)

target_link_libraries(utopia_core utopia_utils utopia_ui)
set_utopia_compiler_features(utopia_core)

install(TARGETS utopia_core
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_core" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")
