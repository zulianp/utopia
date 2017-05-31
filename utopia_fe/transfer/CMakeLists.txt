function(convert_filenames_to_full_paths NAMES)
  unset(tmp_names)
  foreach(name ${${NAMES}})
    list(APPEND tmp_names ${CMAKE_CURRENT_SOURCE_DIR}/${name})
  endforeach()
  set(${NAMES} ${tmp_names} PARENT_SCOPE)
endfunction()



file(GLOB_RECURSE MOONOLITH_SOURCES "*.cpp")
file(GLOB_RECURSE MOONOLITH_HEADERS "*.hpp")
file(GLOB_RECURSE MOONOLITH_HEADERS_2 "*.h")
list(APPEND MOONOLITH_HEADERS ${MOONOLITH_HEADERS_2})

list(APPEND MOONOLITH_INCLUDES 
            moonolith/headers 
            moonolith/external 
            moonolith/kernels 
            ../
            )

convert_filenames_to_full_paths(MOONOLITH_INCLUDES)
list(APPEND MOONOLITH_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})

add_library(moonolith STATIC ${MOONOLITH_SOURCES} ${MOONOLITH_HEADERS})
target_include_directories(moonolith PUBLIC ${MPI_CXX_INCLUDE_PATH} ${MOONOLITH_INCLUDES})
target_link_libraries(moonolith ${TRANSFER_LIB})
target_compile_definitions(moonolith PUBLIC
  "-DWITH_AUTOMATIC_DIFFERENTIATION"
  "-DWITH_MPI"
  "-DPLAIN_INSTALLATION")

list(APPEND HDRS MortarAssembler.hpp)
convert_filenames_to_full_paths(HDRS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(TRANSFER_LIB ${TRANSFER_LIB} moonolith PARENT_SCOPE)
# set(SOURCES )
