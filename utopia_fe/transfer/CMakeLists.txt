list(APPEND TRANSFER_MODULES
  .
  assemblers
  contact_assemblers
  quadrature
)

find_project_files(CMAKE_CURRENT_SOURCE_DIR TRANSFER_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(transfer STATIC ${LOCAL_SOURCES})
target_link_libraries(transfer utopia_fe_core utopia_fe_libmesh_backend utopia_fe_intersection)
target_link_libraries(transfer utopia_fe_visualization)
# target_include_directories(transfer BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

foreach(MODULE ${TRANSFER_MODULES})
	target_include_directories(transfer BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)

set_utopia_compiler_features(transfer)

install(TARGETS transfer
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};transfer" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${MOONOLITH_DEV_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
