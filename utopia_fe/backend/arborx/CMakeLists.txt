# if(UTOPIA_ENABLE_ARBORX)
#     # find_package(MPI COMPONENTS CXX REQUIRED)

#     # FIXME
#     set(MPI_C_FOUND TRUE)
#     set(MPI_C_WORKS TRUE)
#     set(MPI_CXX_WORKS TRUE)
#     set(MPI_CXX_FOUND TRUE)

#     find_package(ArborX REQUIRED)
#     utopiafe_add_module(utopia_arborx ${CMAKE_CURRENT_SOURCE_DIR} ".")
#     target_link_libraries(utopia_arborx ArborX::ArborX Utopia::utopia)


#     set(UTOPIA_FE_LIBRARIES
#         ${UTOPIA_FE_LIBRARIES}
#         PARENT_SCOPE)

# endif()

list(APPEND ARBORX_MODULES .)

if(ArborX_FOUND)
  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${ARBORX_MODULES}"
                   UTOPIA_FE_BUILD_INCLUDES UTOPIA_FE_HEADERS UTOPIA_FE_SOURCES)

  set(UTOPIA_FE_BUILD_INCLUDES
      ${UTOPIA_FE_BUILD_INCLUDES}
      PARENT_SCOPE)

  set(UTOPIA_FE_HEADERS
      ${UTOPIA_FE_HEADERS}
      PARENT_SCOPE)

  set(UTOPIA_FE_SOURCES
      ${UTOPIA_FE_SOURCES}
      PARENT_SCOPE)

endif()