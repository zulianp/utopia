list(APPEND HOMEMADE_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR HOMEMADE_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_fe_homemade STATIC utopia_fe_homemade.cpp ${LOCAL_SOURCES})
target_include_directories(utopia_fe_homemade PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})
foreach(MODULE ${HOMEMADE_MODULES})
	target_include_directories(utopia_fe_homemade PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)

target_link_libraries(utopia_fe_homemade utopia_fe_intersection utopia_fe_core)

install(TARGETS utopia_fe_homemade
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_fe_homemade" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS}")
