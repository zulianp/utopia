list(APPEND INTEROPERABILITY_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR INTEROPERABILITY_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_library(utopia_interoperability STATIC ${LOCAL_SOURCES})

target_link_libraries(utopia_interoperability utopia_core utopia_utils utopia_solvers)

target_include_directories(utopia_interoperability PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(utopia_interoperability PUBLIC ..)
target_include_directories(utopia_interoperability PUBLIC ../backend ${INTEROPERABILITY_MODULES})

# if(TARGET utopia_opencl)
# 	target_link_libraries(utopia_interoperability utopia_opencl)
# endif()

if(TARGET utopia_petsc)
	target_link_libraries(utopia_interoperability utopia_petsc)
endif()

# if(TARGET utopia_blas)
# 	target_link_libraries(utopia_interoperability utopia_blas)
# endif()

# if(TARGET utopia_cuda_cxx)
# 	target_link_libraries(utopia_interoperability utopia_cuda_cxx)
# endif()

if(TARGET utopia_trilinos)
	target_link_libraries(utopia_interoperability utopia_trilinos)
endif()

# if(TARGET utopia_m3elinsol)
# 	target_link_libraries(utopia_interoperability utopia_m3elinsol)
# endif()

install(TARGETS utopia_interoperability
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_interoperability" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")
