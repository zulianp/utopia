
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(MARS)
if(MARS_FOUND)

	list(APPEND MARS_BACKEND_MODULES
		.
	)

	find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_BACKEND_MODULES LOCAL_HEADERS LOCAL_SOURCES)

	add_library(utopia_fe_mars STATIC ${LOCAL_SOURCES} ${LOCAL_HEADERS})
	target_link_libraries(utopia_fe_mars utopia_fe_core utopia_fe_intersection utopia_fe_base)
	target_include_directories(utopia_fe_mars BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	set_utopia_compiler_features(utopia_fe_mars)


	target_include_directories(utopia_fe_mars PUBLIC ${MARS_INCLUDES})
	target_link_libraries(utopia_fe_mars ${MARS_LIBRARIES})

	# set(MARS_INCLUDE_DIR ${MARS_INCLUDE_DIR} PARENT_SCOPE)
	# set(MARS_LIBRARIES ${MARS_LIBRARIES} PARENT_SCOPE)

	#----------------------------------------------------------------------------------------
	#----------------------------------------------------------------------------------------

	install(TARGETS utopia_fe_mars
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		)

	install(FILES ${LOCAL_HEADERS} DESTINATION include)

	set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_mars" PARENT_SCOPE)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

	set(WITH_MARS TRUE PARENT_SCOPE)

else()
	message(WARNING "could not find MARS")
endif()



