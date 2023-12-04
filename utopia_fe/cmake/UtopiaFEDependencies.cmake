# #################UTOPIA###################

find_package(Utopia REQUIRED)
if(Utopia_FOUND)
  message(STATUS "Utopia Found.")
  add_definitions(${UTOPIA_DEFS})

  list(APPEND UTOPIA_FE_BUILD_INCLUDES ${UTOPIA_INCLUDES})
  list(APPEND UTOPIA_FE_DEP_LIBRARIES ${UTOPIA_LIBRARIES})
  list(APPEND UTOPIA_FE_DEP_INCLUDES ${UTOPIA_INCLUDES})


  # Try to get all the variable information from UTOPIA,
  # no sense trying to find them again.
endif()



if(UTOPIA_INSTALL_MOONOLITH)
  include(InstallMoonolith)
  # add_dependencies(utopia_fe par_moonolith)
endif()


# if(UTOPIA_ENABLE_ARBORX)
#   set(ARBORX_SEARCH_PATHS "")

#   if(UTOPIA_ENABLE_ENV_READ)
#     set(ARBORX_SEARCH_PATHS "${ARBORX_SEARCH_PATHS};$ENV{ARBORX_DIR}")
#   endif()

#   find_package(ArborX HINTS ${ARBORX_SEARCH_PATHS} REQUIRED)

#   if(ArborX_FOUND)
#     # Add includes to build_includes, ....
#   endif()
# endif()

if(UTOPIA_ENABLE_LIBMESH)
  find_package(LIBMESH REQUIRED)
  if(LIBMESH_FOUND)
    message(STATUS "Libmesh found.")
    list(APPEND UTOPIA_FE_DEP_LIBRARIES ${LIBMESH_LIBRARIES})
    list(APPEND UTOPIA_FE_DEP_INCLUDES ${LIBMESH_INCLUDE_DIR})
  endif()
endif()

if(UTOPIA_ENABLE_MOONOLITH)
  find_package(ParMoonolith QUIET)
  if(ParMoonolith_FOUND)
    message(STATUS "ParMoonolith found.")
    list(APPEND UTOPIA_FE_DEP_LIBRARIES ${MOONOLITH_LIBRARIES})
    list(APPEND UTOPIA_FE_DEP_INCLUDES ${MOONOLITH_INCLUDES})
  endif()
endif()

