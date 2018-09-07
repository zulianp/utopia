
# find dependencies
find_package(Trilinos)

if(TRILINOS_FOUND)
    foreach(LIB ${Trilinos_LIBRARY_DIRS})
        list(APPEND UTOPIA_LIBRARIES -L${LIB})
    endforeach(LIB)

    #	list(APPEND UTOPIA_INCLUDES ${Trilinos_INCLUDE_DIRS})
    # list(APPEND UTOPIA_LIBRARIES ${Trilinos_LIBRARY_DIRS})
    list(APPEND Trilinos_all_libs ${Trilinos_TPL_LIBRARIES})
    list(APPEND Trilinos_all_libs ${Trilinos_LIBRARIES})
    list(REMOVE_DUPLICATES Trilinos_all_libs)
    foreach(LIB ${Trilinos_all_libs})
        if(EXISTS ${LIB})
            list(APPEND UTOPIA_LIBRARIES ${LIB})
        else()
            list(APPEND UTOPIA_LIBRARIES -l${LIB})
        endif()
    endforeach()
    list(APPEND Trilinos_TPL_INCLUDE_DIRS ${Trilinos_additional_headers})
    list(APPEND UTOPIA_INCLUDES ${Trilinos_TPL_INCLUDE_DIRS})
    #####
    message(STATUS "Trilinos_additional_headers ${Trilinos_additional_headers}")
    message(STATUS "Trilinos_TPL_INCLUDE_DIRS ${Trilinos_TPL_INCLUDE_DIRS})

    # set(CMAKE_C_COMPILER ${TRILINOS_COMPILER})

    if(NOT MPI_CXX_COMPILER)
		set(MPI_CXX_COMPILER $ENV{MPI_CXX_COMPILER})
		message(STATUS "compiler ${MPI_CXX_COMPILER}")
    endif()

	if(MPI_CXX_COMPILER)
		set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
		set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
	else()
		 execute_process(COMMAND mpicxx -v RESULT_VARIABLE MPICXX_FAILED)

		if(MPICXX_FAILED)
			message(STATUS "Using CMAKE compiler, you can define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>")
		else()
    		message(STATUS "-----------------------------------------------")
	    	message(STATUS "\n[MPI] using mpicxx for compiling c++ files.\nIf you want to use your own compiler define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>")
	    	message(STATUS "-----------------------------------------------")
	    	set(CMAKE_CXX_COMPILER mpicxx)
	    	set(CMAKE_CXX_COMPILER_DEBUG mpicxx)
		endif()
	endif()

    set(WITH_TRILINOS TRUE)
    if(KOKKOS_CUDA_HPP)
        set(KOKKOS_CUDA TRUE)
    endif()
    if(KOKKOS_OPENMP_HPP)
        set(KOKKOS_OPENMP TRUE)
    endif()

else()
    message(WARNING "[Warning] Trilinos not found")
    set(WITH_TRILINOS FALSE)
endif()

# set-up module
if(TRILINOS_FOUND)
	set(UTOPIA_TRILINOS TRUE PARENT_SCOPE)
	set(WITH_TRILINOS ON PARENT_SCOPE)

	list(APPEND TRILINOS_MODULES
		.
		solvers
		solvers/linear
		preconditioners
		)

	find_project_files(CMAKE_CURRENT_SOURCE_DIR TRILINOS_MODULES LOCAL_HEADERS LOCAL_SOURCES)
	add_library(utopia_trilinos STATIC utopia_trilinos.cpp ${LOCAL_SOURCES})

	foreach(MODULE ${TRILINOS_MODULES})
		target_include_directories(utopia_trilinos PUBLIC ${MODULE})
	endforeach(MODULE)

	target_include_directories(utopia_trilinos PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}  ${CMAKE_CURRENT_SOURCE_DIR}/..)
	target_include_directories(utopia_trilinos PUBLIC ${Trilinos_INCLUDE_DIRS})

	target_link_libraries(utopia_trilinos utopia_core utopia_utils utopia_solvers)
	target_link_libraries(utopia_trilinos ${Trilinos_LIBRARIES})
	target_link_libraries(utopia_trilinos ${Trilinos_TPL_LIBRARIES})
	target_link_libraries(utopia_trilinos "-L${Trilinos_LIBRARY_DIRS}")
	add_definitions(${TRILINOS_DEFINITIONS})


	install(TARGETS utopia_trilinos
		RUNTIME DESTINATION bin 
		LIBRARY DESTINATION lib 
		ARCHIVE DESTINATION lib 
		)

	install(FILES ${LOCAL_HEADERS} DESTINATION include)

	set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_trilinos" PARENT_SCOPE)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")


	set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
	set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
	set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
endif()
