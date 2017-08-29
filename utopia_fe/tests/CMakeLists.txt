

list(APPEND LOCAL_HEADERS 
	utopia_SemigeometricMultigridTest.hpp
	utopia_UtopiaFETests.hpp
	)

list(APPEND LOCAL_SOURCES 
	utopia_SemigeometricMultigridTest.cpp
	utopia_UtopiaFETests.cpp
	)

add_library(utopia_fe_tests STATIC ${LOCAL_SOURCES})
target_link_libraries(utopia_fe_tests transfer)
target_include_directories(utopia_fe_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) #;${CMAKE_CURRENT_SOURCE_DIR}/../transfer;${CMAKE_CURRENT_SOURCE_DIR}/.. )


install(TARGETS utopia_fe_tests
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_tests" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -g")
