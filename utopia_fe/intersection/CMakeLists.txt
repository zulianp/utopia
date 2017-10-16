
set(LOCAL_SOURCES 
  utopia_triangulate.cpp
  utopia_intersector.cpp
  clipper/clipper.cpp
  )

set(LOCAL_HEADERS
  utopia_triangulate.hpp
  utopia_intersector.hpp
  clipper/clipper.hpp
)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

find_package(ParMoonolith)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

# include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(utopia_fe_intersection STATIC ${LOCAL_SOURCES})
add_dependencies(utopia_fe_intersection par_moonolith)
target_include_directories(utopia_fe_intersection PUBLIC ${MOONOLITH_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/clipper)
target_link_libraries(utopia_fe_intersection ${MOONOLITH_LIB})
target_link_libraries(utopia_fe_intersection ${MOONOLITH_LIB})

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};utopia_fe_intersection" PARENT_SCOPE)

set(MOONOLITH_LIB ${MOONOLITH_LIB} PARENT_SCOPE)
set(MOONOLITH_INCLUDES ${MOONOLITH_INCLUDES} PARENT_SCOPE)

install(TARGETS utopia_fe_intersection
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${LOCAL_HEADERS} DESTINATION include)