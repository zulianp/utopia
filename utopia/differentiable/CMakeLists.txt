list(APPEND DIFFERENTIABLE_MODULES
	.
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR DIFFERENTIABLE_MODULES LOCAL_HEADERS LOCAL_SOURCES)

list(APPEND LOCAL_SOURCES utopia_differentiable.cpp)

add_library(utopia_differentiable STATIC ${LOCAL_SOURCES})
target_link_libraries(utopia_differentiable utopia_core)
target_include_directories(utopia_differentiable PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) 


install(TARGETS utopia_differentiable
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_differentiable" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -g")
