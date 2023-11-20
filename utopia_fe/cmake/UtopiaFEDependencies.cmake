# #################UTOPIA###################

set(UTOPIA_SEARCH_PATHS "/usr/local/;/usr/")

if(UTOPIA_ENABLE_ENV_READ)
  set(UTOPIA_SEARCH_PATHS "${UTOPIA_SEARCH_PATHS};$ENV{UTOPIA_DIR};$ENV{Utopia_DIR}")
endif()

find_package(Utopia HINTS ${UTOPIA_SEARCH_PATHS} REQUIRED)
if(Utopia_FOUND)
  message(STATUS "Utopia Found.")
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


