list(APPEND TEST_MODULES
	.
	simulations
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR TEST_MODULES LOCAL_HEADERS LOCAL_SOURCES)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_core" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")


add_library(utopia_fe_tests STATIC ${LOCAL_SOURCES})
target_link_libraries(
	utopia_fe_tests
	transfer
	utopia_fe_visualization
	utopia_fe_homemade
	utopia_fe_libmesh_backend
	utopia_fe_contact
	utopia_fe_materials
	utopia_fe_multigrid
	utopia_fe_functions
	utopia_fe_fractures
)

target_include_directories(utopia_fe_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/simulations) 
set_utopia_compiler_features(utopia_fe_tests)

install(TARGETS utopia_fe_tests
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_tests" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
