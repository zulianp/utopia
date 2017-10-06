
list(APPEND LOCAL_HEADERS 
	utopia_Backends.hpp
	utopia_ScalarBackend.hpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

install(FILES ${LOCAL_HEADERS} DESTINATION include)
