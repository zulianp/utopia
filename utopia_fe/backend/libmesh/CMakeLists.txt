list(APPEND LIBMESH_BACKEND_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR LIBMESH_BACKEND_MODULES LOCAL_HEADERS LOCAL_SOURCES)


add_library(utopia_fe_libmesh_backend STATIC ${LOCAL_SOURCES})
target_link_libraries(utopia_fe_libmesh_backend utopia_fe_core utopia_fe_intersection utopia_fe_base)
target_include_directories(utopia_fe_libmesh_backend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_utopia_compiler_features(utopia_fe_libmesh_backend)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(LIBMESH)

if(LIBMESH_FOUND)
	target_include_directories(utopia_fe_libmesh_backend PUBLIC ${LIBMESH_INCLUDE_DIRS})
	target_link_libraries(utopia_fe_libmesh_backend ${LIBMESH_LIBRARIES})

	set(LIBMESH_INCLUDE_DIR ${LIBMESH_INCLUDE_DIR} PARENT_SCOPE)
	set(LIBMESH_LIBRARIES ${LIBMESH_LIBRARIES} PARENT_SCOPE)
else()
	message(FATAL_ERROR "could not find LibMesh")
endif()

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

install(TARGETS utopia_fe_libmesh_backend
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_libmesh_backend" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
