
list(APPEND LOCAL_HEADERS 
	utopia_LibMeshBackend.hpp
	utopia_LibMeshContext.hpp
	utopia_LibMeshLambdaAssembly.hpp
	utopia_NormalTangentialCoordinateSystem.hpp
	utopia_UGMeshReader.hpp
	)

list(APPEND LOCAL_SOURCES 
	utopia_LibMeshBackend.cpp
	utopia_NormalTangentialCoordinateSystem.cpp
	utopia_UGMeshReader.cpp
	)

add_library(utopia_fe_libmesh_backend STATIC ${LOCAL_SOURCES})
target_link_libraries(utopia_fe_libmesh_backend utopia_fe_core utopia_fe_intersection)
target_include_directories(utopia_fe_libmesh_backend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(LIBMESH)

if(LIBMESH_FOUND)
	target_include_directories(utopia_fe_libmesh_backend PUBLIC ${LIBMESH_INCLUDE_DIR})
	target_link_libraries(utopia_fe_libmesh_backend ${LIBMESH_LIBRARIES})
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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${MOONOLITH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
