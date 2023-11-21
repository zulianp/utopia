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

# #################MPI######################
if(UTOPIA_ENABLE_MPI)

  find_package(MPIExtended REQUIRED)
  if(MPI_FOUND)
    if(MPI_C_INCLUDE_PATH)
      set(UTOPIA_DEP_INCLUDES "${UTOPIA_DEP_INCLUDES};${MPI_C_INCLUDE_PATH}")
    endif()

    if(MPI_CXX_INCLUDE_PATH)
      set(UTOPIA_DEP_INCLUDES "${UTOPIA_DEP_INCLUDES};${MPI_CXX_INCLUDE_PATH}")
    endif()

    if(MPI_LIBRARIES)
      set(UTOPIA_DEP_LIBRARIES "${UTOPIA_DEP_LIBRARIES};${MPI_LIBRARIES}")
    endif()

    if(MPI_C_LIBRARIES)
      set(UTOPIA_DEP_LIBRARIES "${UTOPIA_DEP_LIBRARIES};${MPI_C_LIBRARIES}")
    endif()

    if(MPI_CXX_LIBRARIES)
      set(UTOPIA_DEP_LIBRARIES "${UTOPIA_DEP_LIBRARIES};${MPI_CXX_LIBRARIES}")

      set(UTOPIA_MPI_DIR ${MPI_LIBRARIES})
      set(UTOPIA_MPI_VERSION ${MPI_C_VERSION})
    endif()
  else()
    message(WARNING "NO Proper MPI installation")
  endif()
endif()

if(MPI_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
endif()

if(UTOPIA_ENABLE_ARBORX)
  set(ARBORX_SEARCH_PATHS "")

  if(UTOPIA_ENABLE_ENV_READ)
    set(ARBORX_SEARCH_PATHS "${ARBORX_SEARCH_PATHS};$ENV{ARBORX_DIR}")
  endif()

  find_package(ArborX HINTS ${ARBORX_SEARCH_PATHS} REQUIRED)

  if(ArborX_FOUND)
    # Add includes to build_includes, ....
  endif()
endif()

if(UTOPIA_ENABLE_LIBMESH)
  include(InstallLIBMESH)
  find_package(LIBMESH REQUIRED)

  if(LIBMESH_FOUND)
    message(STATUS "Libmesh found.")
  endif()
endif()

if(UTOPIA_ENABLE_MOONOLITH)
  find_package(ParMoonolith REQUIRED)
  if(ParMoonolith_FOUND)
    message(STATUS "ParMoonolith found.")
  endif()
endif()

if(UTOPIA_ENABLE_TRILINOS)
  # set(TRILINOS_SEARCH_PATHS "/usr/;/usr/local/")

  # if(UTOPIA_ENABLE_ENV_READ)
  #   set(TRILINOS_SEARCH_PATHS
  #       "${TRILINOS_SEARCH_PATHS};$ENV{TRILINOS_DIR};$ENV{TRILINOS_DIR}/lib/cmake/Trilinos"
  #   )
  # endif()

  find_package(
    Trilinos
    NO_MODULE
    NO_CMAKE_ENVIRONMENT_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    HINTS
    ${TRILINOS_SEARCH_PATHS}
    REQUIRED)

  if(Trilinos_FOUND)

    foreach(LIB ${Trilinos_LIBRARY_DIRS})
      list(APPEND UTOPIA_THIRDPARTY_LIBRARIES -L${LIB})
    endforeach(LIB)

    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${Trilinos_INCLUDE_DIRS})
    # list(APPEND UTOPIA_LIBRARIES ${Trilinos_LIBRARY_DIRS})
    list(APPEND Trilinos_all_libs ${Trilinos_LIBRARIES})
    list(APPEND Trilinos_all_libs ${Trilinos_TPL_LIBRARIES})
    list(REVERSE Trilinos_all_libs)
    list(REMOVE_DUPLICATES Trilinos_all_libs)
    list(REVERSE Trilinos_all_libs)
    foreach(LIB ${Trilinos_all_libs})
      if(EXISTS ${LIB})
        list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${LIB})
      else()
        list(APPEND UTOPIA_THIRDPARTY_LIBRARIES -l${LIB})
      endif()
    endforeach()
    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${Trilinos_TPL_INCLUDE_DIRS})
    # ############## End of population of variables for utopia-config.makefile
    # #########################3

    set(UTOPIA_TRILINOS_DEPS Kokkos::kokkos)
    # Unfortunately trilinos libraries are not namespaced (yet), we use the bare
    # target names
    list(APPEND Utopia_Trilinos_possible_packages Amesos2 Belos Ifpack2 MueLu)
    foreach(package ${Utopia_Trilinos_possible_packages})
      # e.g Amesos2: UTOPIA_ENABLE_TRILINOS_AMESOS2=TRUE if package is found in
      # Trilinos_PACKAGE_LIST
      string(TOUPPER ${package} packageUpper)
      string(TOLOWER ${package} packageLower)
      list(FIND Trilinos_PACKAGE_LIST ${package} PACKAGE_FOUND)
      if(NOT PACKAGE_FOUND EQUAL -1)
        set(UTOPIA_ENABLE_TRILINOS_${packageUpper} TRUE)
        if(${package} STREQUAL "MueLu")
          # ugly hack, but we need to link with muelu-adapters also I cannot #
          # wait until Trilinos finally supports cmake targets correctly
          list(APPEND UTOPIA_TRILINOS_DEPS ${MueLu_LIBRARIES})
        endif()
        # message(STATUS "${package}")
        list(APPEND UTOPIA_TRILINOS_DEPS ${${package}_LIBRARIES})
      endif()
    endforeach()

    # CHECK UTOPIA_THIRDPARTY_LIBRARIES: it is used in the makefile lines in
    # main cmake.
    list(APPEND UTOPIA_DEP_LIBRARIES ${UTOPIA_TRILINOS_DEPS})

    find_package(TpetraExt)
    if(TRILINOS_TPETRAEXT_FOUND)
      set(UTOPIA_ENABLE_TRILINOS_TPETRAEXT TRUE)
    endif()

    set(UTOPIA_TRILINOS_DIR ${Trilinos_DIR})
    set(UTOPIA_TRILINOS_VERSION ${Trilinos_VERSION})

  else()
    message(WARNING "[Warning] Trilinos not found")
  endif()
  add_subdirectory(backend/trilinos)

endif()
