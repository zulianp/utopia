list(APPEND CORE_MODULES
	.
	eval
	expressions
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR CORE_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_core STATIC ${LOCAL_SOURCES})
target_include_directories(utopia_fe_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_utopia_compiler_features(utopia_fe_core)

foreach(MODULE ${CORE_MODULES})
	target_include_directories(utopia_fe_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)


install(TARGETS utopia_fe_core
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_core" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
