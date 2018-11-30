
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(Intrepid2)
if(INTREPID2_FOUND)

	list(APPEND INTREPID2_BACKEND_MODULES
		.
	)

	find_project_files(CMAKE_CURRENT_SOURCE_DIR INTREPID2_BACKEND_MODULES LOCAL_HEADERS LOCAL_SOURCES)


	add_library(utopia_fe_intrepid2 STATIC ${LOCAL_SOURCES} ${LOCAL_HEADERS})
	target_link_libraries(utopia_fe_intrepid2 utopia_fe_core utopia_fe_intersection utopia_fe_base)
	target_include_directories(utopia_fe_intrepid2 BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	set_utopia_compiler_features(utopia_fe_intrepid2)


	target_include_directories(utopia_fe_intrepid2 PUBLIC ${INTREPID2_INCLUDE_DIR})
	target_link_libraries(utopia_fe_intrepid2 utopia_fe_libmesh_backend ${INTREPID2_LIBRARIES})

	# set(INTREPID2_INCLUDE_DIR ${Intrepid2_INCLUDE_DIR} PARENT_SCOPE)
	# set(Intrepid2_LIBRARIES ${Intrepid2_LIBRARIES} PARENT_SCOPE)

	#----------------------------------------------------------------------------------------
	#----------------------------------------------------------------------------------------

	install(TARGETS utopia_fe_intrepid2
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		)

	install(FILES ${LOCAL_HEADERS} DESTINATION include)

	set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_intrepid2" PARENT_SCOPE)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

	set(WITH_INTREPID2 TRUE PARENT_SCOPE)

else()
	message(WARNING "could not find Intrepid2")
endif()



