
# find dependency
find_package(CUDA)
if(CUDA_FOUND)
	find_package(Thrust)
	if(THRUST_FOUND)
		file(GLOB CUDA_FILES "*.cu") 
		set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) 
		cuda_include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${THRUST_INCLUDE_DIR} ${UTOPIA_ROOT_PATH}/core ${UTOPIA_ROOT_PATH}/core/expressions ${CMAKE_BINARY_DIR})
		set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; "-std=c++11")
		cuda_compile(CUDA_OBJECTS ${CUDA_FILES})
		message(STATUS "${CUDA_OBJECTS}")
		cuda_add_library(utopia_cuda ${CUDA_FILES})
		set(WITH_CUDA ON)
	else()
		message(FATAL_ERROR "Thrust not found")	
	endif()
else()
	message(FATAL_ERROR "Cuda not found")
endif()


# set-up module
if(THRUST_FOUND)
	set(UTOPIA_CUDA TRUE PARENT_SCOPE)
	set(WITH_CUDA ON PARENT_SCOPE)

	list(APPEND CUDA_MODULES
		.
		)

	find_project_files(CMAKE_CURRENT_SOURCE_DIR CUDA_MODULES LOCAL_HEADERS LOCAL_SOURCES)

	add_library(utopia_cuda_cxx STATIC utopia_cuda.cpp ${LOCAL_SOURCES})

	foreach(MODULE ${CUDA_MODULES})
		target_include_directories(utopia_cuda_cxx PUBLIC ${MODULE})
	endforeach(MODULE)

	target_include_directories(utopia_cuda_cxx PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}  ${CMAKE_CURRENT_SOURCE_DIR}/..)
	target_include_directories(utopia_cuda_cxx PUBLIC ${CUDA_INCLUDE_DIRS})

	target_link_libraries(utopia_cuda_cxx utopia_core utopia_utils utopia_solvers utopia_cuda)
	target_link_libraries(utopia_cuda_cxx ${CUDA_LIBRARIES})

	install(TARGETS utopia_cuda_cxx
		RUNTIME DESTINATION bin 
		LIBRARY DESTINATION lib 
		ARCHIVE DESTINATION lib 
		)

	install(TARGETS utopia_cuda
		RUNTIME DESTINATION bin 
		LIBRARY DESTINATION lib 
		ARCHIVE DESTINATION lib 
		)

	install(FILES ${LOCAL_HEADERS} DESTINATION include)

	set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_cuda_cxx;utopia_cuda" PARENT_SCOPE)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")

	set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
	set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
	set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
endif()
