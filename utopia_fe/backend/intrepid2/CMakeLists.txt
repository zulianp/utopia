if(Intrepid2_FOUND)
    list(APPEND INTREPID_MODULES .)

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${INTREPID_MODULES}"
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


# find_package(Intrepid2 HINTS $ENV{TRILINOS_DIR}/lib/cmake/Intrepid2 ${Trilinos_DIR}/../Intrepid2 REQUIRED)

# utopiafe_add_module(utopia_intrepid2 ${CMAKE_CURRENT_SOURCE_DIR} ".")

# target_link_directories(utopia_intrepid2 PUBLIC ${Intrepid2_LIBRARY_DIRS})
# target_link_libraries(utopia_intrepid2 PUBLIC utopia_fe_core utopia_fe_frontend utopia_fe_kokkos Utopia::utopia)
# target_link_libraries(utopia_intrepid2 PUBLIC ${Intrepid2_LIBRARIES})

# target_include_directories(utopia_intrepid2
#                            PUBLIC ${Intrepid2_TPL_INCLUDE_DIRS})

# target_include_directories(utopia_intrepid2 PUBLIC ${Intrepid2_INCLUDE_DIRS})

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DHAVE_INTREPID2_DEBUG")
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" PARENT_SCOPE)

# set(UTOPIA_ENABLE_INTREPID2 TRUE)

# set(UTOPIA_ENABLE_INTREPID2
#     TRUE
#     PARENT_SCOPE)

# set(UTOPIA_FE_LIBRARIES
#     ${UTOPIA_FE_LIBRARIES}
#     PARENT_SCOPE)
