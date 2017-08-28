function(convert_filenames_to_full_paths NAMES)
  unset(tmp_names)
  foreach(name ${${NAMES}})
    list(APPEND tmp_names ${CMAKE_CURRENT_SOURCE_DIR}/${name})
  endforeach()
  set(${NAMES} ${tmp_names} PARENT_SCOPE)
endfunction()

find_package(ParMoonolith)
message(STATUS "Ml: ${MOONOLITH_LIB}")

set(TRANSFER_SOURCES 
  utopia_BoxAdapter.cpp
  utopia_assemble_contact.cpp
  Box.cpp
  HashGrid.cpp
  LibmeshMultigridForMoose.cpp
  MortarAssemble.cpp
  MortarAssembler.cpp
  ParMortarAssembler.cpp
  opencl_adapter.cpp
  triangulate.cpp
  utopia_BoxAdapter.cpp
  utopia_FESpaceAdapter.cpp
  utopia_FESpacesAdapter.cpp
  utopia_FESpacesRAdapter.cpp
  utopia_Socket.cpp
  utopia_assemble_contact.cpp
  utopia_assemble_volume_transfer.cpp
  utopia_assemble_volume_transfer_r.cpp
  utopia_copy_dofmap.cpp
  )

list(APPEND MOONOLITH_INCLUDES 
            ${CMAKE_CURRENT_SOURCE_DIR}
            ../ )

add_library(transfer STATIC ${TRANSFER_SOURCES})
add_dependencies(transfer par_moonolith)
target_include_directories(transfer PUBLIC ${MOONOLITH_INCLUDES})
target_link_libraries(transfer ${MOONOLITH_LIB})


list(APPEND HDRS MortarAssembler.hpp utopia_assemble_contact.hpp)

convert_filenames_to_full_paths(HDRS)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(TRANSFER_LIB ${TRANSFER_LIB} transfer PARENT_SCOPE)

