list(APPEND CORE_MODULES
	.
	eval
	expressions
	expressions_parallel
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR CORE_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_core STATIC utopia_core.cpp ${LOCAL_SOURCES})
target_include_directories(utopia_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})
foreach(MODULE ${CORE_MODULES})
	target_include_directories(utopia_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)

target_link_libraries(utopia_core utopia_utils)

install(TARGETS utopia_core
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_LIBRARIES "${UTOPIA_LIBRARIES};utopia_core" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS}")
