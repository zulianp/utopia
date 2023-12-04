# #################UTOPIA###################

find_package(Utopia REQUIRED)
if(Utopia_FOUND)
  message(STATUS "Utopia Found.")
  add_definitions(${UTOPIA_DEFS})

  list(APPEND UTOPIA_FE_BUILD_INCLUDES ${UTOPIA_INCLUDES})
  list(APPEND UTOPIA_FE_DEP_LIBRARIES ${UTOPIA_LIBRARIES})
  list(APPEND UTOPIA_FE_DEP_INCLUDES ${UTOPIA_INCLUDES})

  # Try to get all the variable information from UTOPIA, no sense trying to find
  # them again.
endif()

if(UTOPIA_INSTALL_MOONOLITH)
  include(InstallMoonolith)
  # add_dependencies(utopia_fe par_moonolith)
endif()

# if(UTOPIA_ENABLE_ARBORX) set(ARBORX_SEARCH_PATHS "")

# if(UTOPIA_ENABLE_ENV_READ) set(ARBORX_SEARCH_PATHS
# "${ARBORX_SEARCH_PATHS};$ENV{ARBORX_DIR}") endif()

# find_package(ArborX HINTS ${ARBORX_SEARCH_PATHS} REQUIRED)

# if(ArborX_FOUND) # Add includes to build_includes, .... endif() endif()

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

if(UTOPIA_ENABLE_INTREPID2)
  set(UTOPIA_ENABLE_INTREPID2_SEARCH_PATHS "${Trilinos_DIR}/../Intrepid2")
  if(UTOPIA_ENABLE_ENV_READ)
    set(UTOPIA_ENABLE_INTREPID2_SEARCH_PATHS
        "${UTOPIA_ENABLE_INTREPID2_SEARCH_PATHS};$ENV{TRILINOS_DIR}/lib/cmake/Intrepid2"
    )
  endif()
  find_package(Intrepid2 HINTS ${UTOPIA_ENABLE_INTREPID2_SEARCH_PATHS} REQUIRED)
  if(Intrepid2_FOUND)
    message(STATUS "Intrepid2 found.")
    # list(APPEND UTOPIA_FE_DEP_LIBRARIES ${MOONOLITH_LIBRARIES}) list(APPEND
    # UTOPIA_FE_DEP_INCLUDES ${MOONOLITH_INCLUDES})
  endif()

  find_package(
    Kokkos
    NO_MODULE
    NO_CMAKE_ENVIRONMENT_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    HINTS
    $ENV{KOKKOS_DIR}
    $ENV{TRILINOS_DIR}/lib/cmake/Kokkos
    ${Trilinos_DIR}/../Kokkos
    REQUIRED)
endif()

if(UTOPIA_ENABLE_MARS)
  find_package(Mars 0 REQUIRED)

  message(STATUS "MARS_LIBRARIES=${Mars_LIBRARIES}")
  if(UTOPIA_ENABLE_MARS_VTK)
    find_package(
      VTK
      COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersGeneral vtkIOXML
                 vtkIOParallel vtkIOParallelXML
      REQUIRED)
  endif()
endif()
