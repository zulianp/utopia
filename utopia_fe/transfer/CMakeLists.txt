
list(APPEND LOCAL_HEADERS 
  Box.hpp
  HashGrid.hpp
  LibmeshContactForMoose.hpp
  LibmeshMultigridForMoose.hpp
  LibmeshTransferForMoose.hpp
  LibmeshTransferForMooseReverse.hpp
  MortarAssemble.hpp
  MortarAssembler.hpp
  ParMortarAssembler.hpp
  triangulate.hpp
  utopia_BoxAdapter.hpp
  utopia_ElementDofMap.hpp
  utopia_FESpaceAdapter.hpp
  utopia_FESpacesAdapter.hpp
  utopia_FESpacesRAdapter.hpp
  utopia_Polygon.hpp
  utopia_SElementAdapter.hpp
  utopia_STree.hpp
  utopia_VElementAdapter.hpp
  utopia_VTree.hpp
  utopia_assemble_contact.hpp
  utopia_assemble_volume_transfer.hpp
  utopia_assemble_volume_transfer_r.hpp
  utopia_copy_dofmap.hpp
  )

list(APPEND LOCAL_SOURCES 
  Box.cpp
  HashGrid.cpp
  LibmeshMultigridForMoose.cpp
  MortarAssemble.cpp
  MortarAssembler.cpp
  ParMortarAssembler.cpp
  utopia_BoxAdapter.cpp
  utopia_FESpaceAdapter.cpp
  utopia_FESpacesAdapter.cpp
  utopia_FESpacesRAdapter.cpp
  utopia_assemble_contact.cpp
  utopia_assemble_volume_transfer.cpp
  utopia_assemble_volume_transfer_r.cpp
  utopia_copy_dofmap.cpp
  )

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
