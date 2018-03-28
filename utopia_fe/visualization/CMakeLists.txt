
list(APPEND LOCAL_HEADERS 
	utopia_Socket.hpp
	)

list(APPEND LOCAL_SOURCES 
	utopia_Socket.cpp
	)

add_library(utopia_fe_visualization STATIC ${LOCAL_SOURCES})
target_link_libraries(utopia_fe_visualization utopia_fe_libmesh_backend transfer)
target_include_directories(utopia_fe_visualization PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})
set_utopia_compiler_features(utopia_fe_visualization)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
if(ENABLE_BOOST)
	find_package(Boost)
	find_package(BoostExtended)
	if(Boost_FOUND)
			target_include_directories(utopia_fe_visualization PRIVATE ${Boost_INCLUDE_DIRS})
			link_directories(${Boost_LIBRARY_DIRS})
			target_link_libraries(utopia_fe_visualization ${Boost_LIBRARIES})
			# add_definitions("-DWITH_BOOST")
			set(WITH_BOOST TRUE PARENT_SCOPE)
	else()
		MESSAGE(WARNING "Boost not found.")
		set(WITH_BOOST FALSE PARENT_SCOPE)
	endif()
else()
	set(WITH_BOOST FALSE PARENT_SCOPE)
endif()
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

install(TARGETS utopia_fe_visualization
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_visualization" PARENT_SCOPE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
