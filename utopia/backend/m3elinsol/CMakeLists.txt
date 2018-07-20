find_package(M3ELINSOL)

if(M3ELINSOL_FOUND)
    list(APPEND UTOPIA_INCLUDES ${M3ELINSOL_INCLUDES})
    list(APPEND UTOPIA_LIBRARIES ${M3ELINSOL_LIBRARIES})
    list(APPEND UTOPIA_DEFS ${M3ELINSOL_DEFS})
    set(WITH_M3ELINSOL ON)
else()
    message(WARNING "[Warning] M3ELINSOL not found")
    set(WITH_M3ELINSOL OFF)
endif()

# set-up module
if(M3ELINSOL_FOUND)
	set(UTOPIA_M3ELINSOL TRUE PARENT_SCOPE)
	set(WITH_M3ELINSOL ON PARENT_SCOPE)

	list(APPEND M3ELINSOL_MODULES
		.
	)

	


	find_project_files(CMAKE_CURRENT_SOURCE_DIR M3ELINSOL_MODULES LOCAL_HEADERS LOCAL_SOURCES)
# compiling also m3elinsol c++ interface file
	list(APPEND LOCAL_SOURCES
			${M3ELINSOL_SOURCE_FILES}
	)
# message(STATUS "${LOCAL_SOURCES}")
	add_library(utopia_m3elinsol STATIC ${LOCAL_SOURCES})

	foreach(MODULE ${M3ELINSOL_MODULES})
		target_include_directories(utopia_m3elinsol PUBLIC ${MODULE})
	endforeach(MODULE)

	target_include_directories(utopia_m3elinsol PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}  ${CMAKE_CURRENT_SOURCE_DIR}/..)
	target_include_directories(utopia_m3elinsol PUBLIC ${M3ELINSOL_INCLUDES})

	target_link_libraries(utopia_m3elinsol utopia_core utopia_utils utopia_solvers)
	target_link_libraries(utopia_m3elinsol ${M3ELINSOL_LIBRARIES})

	find_package(OpenMP)
	if(OPENMP_FOUND)
		SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}" PARENT_SCOPE)
		SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" PARENT_SCOPE)
		SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" PARENT_SCOPE)
		set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${OpenMP_CXX_FLAGS} -lgfortran" PARENT_SCOPE)
	endif()

	find_package(gfortran)
	if(GFORTRAN_FOUND)
		message(STATUS "Found gfortran")
	endif()

	# set_target_properties(utopia_m3elinsol PROPERTIES COMPILE_FLAGS "-fopenmp")
	# set_target_properties(utopia_m3elinsol PROPERTIES LINK_FLAGS "-fopenmp")

	if(TARGET utopia_petsc)
		target_link_libraries(utopia_m3elinsol utopia_petsc)
	endif()

	if(TARGET utopia_blas)
		target_link_libraries(utopia_m3elinsol utopia_blas)
	endif()

	if(TARGET utopia_trilinos)
		target_link_libraries(utopia_m3elinsol utopia_trilinos)
	endif()

	add_definitions(${M3ELINSOL_DEFS})

	install(TARGETS utopia_m3elinsol
		RUNTIME DESTINATION bin 
		LIBRARY DESTINATION lib 
		ARCHIVE DESTINATION lib 
	)

	install(FILES ${LOCAL_HEADERS} DESTINATION include)

	set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_m3elinsol" PARENT_SCOPE)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")

	set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
	set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
	set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
endif()
