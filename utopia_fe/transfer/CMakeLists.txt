list(APPEND TRANSFER_MODULES
  .
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR TRANSFER_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(transfer STATIC ${LOCAL_SOURCES})
target_link_libraries(transfer utopia_fe_core utopia_fe_libmesh_backend utopia_fe_intersection)

target_link_libraries(transfer utopia_fe_visualization) 

target_include_directories(transfer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS transfer
  RUNTIME DESTINATION bin 
  LIBRARY DESTINATION lib 
  ARCHIVE DESTINATION lib 
  )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};transfer" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${MOONOLITH_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
