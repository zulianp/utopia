list(APPEND TEST_MODULES
	.
	test_problems
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR TEST_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_library(utopia_tests STATIC ${LOCAL_SOURCES})

target_link_libraries(utopia_tests utopia_core utopia_utils utopia_solvers utopia_differentiable)
target_include_directories(utopia_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) 
target_include_directories(utopia_tests PUBLIC ..)
target_include_directories(utopia_tests PUBLIC ../backend) 

if(TARGET utopia_opencl)
	target_link_libraries(utopia_tests utopia_opencl)
endif()

if(TARGET utopia_petsc)
	target_link_libraries(utopia_tests utopia_petsc)
endif()

if(TARGET utopia_blas)
	target_link_libraries(utopia_tests utopia_blas)
endif()

if(TARGET utopia_cuda)
	target_link_libraries(utopia_tests utopia_cuda)
endif()

install(TARGETS utopia_tests
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_LIBRARIES "${UTOPIA_LIBRARIES};utopia_tests" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -g")
