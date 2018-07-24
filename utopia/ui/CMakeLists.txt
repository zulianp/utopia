list(APPEND UI_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR UI_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_ui STATIC ${LOCAL_SOURCES})
target_include_directories(utopia_ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})
target_include_directories(utopia_ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(utopia_ui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

#hidden dependencies
target_include_directories(utopia_ui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rapidxml)
target_include_directories(utopia_ui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rapidjson)

if(WIN32)
	target_include_directories(utopia_ui PRIVATE windows)
endif()


target_link_libraries(utopia_ui utopia_utils)

set_utopia_compiler_features(utopia_ui)

install(TARGETS utopia_ui
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_ui" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")


set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
