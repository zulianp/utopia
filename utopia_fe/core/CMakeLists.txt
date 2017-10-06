
list(APPEND LOCAL_HEADERS 
	utopia_MixedFESpace.hpp
	utopia_fe_core.hpp
	utopia_fe_lang.hpp
	utopia_block_var.hpp
	)

list(APPEND LOCAL_SOURCES 
	utopia_fe_core.cpp
	)

add_library(utopia_fe_core STATIC ${LOCAL_SOURCES})
target_include_directories(utopia_fe_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS utopia_fe_core
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_core" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${MOONOLITH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
