# ##############################################################################
# LOCAL INSTALLS

if(UTOPIA_INSTALL_MOONOLITH)
  include(InstallMoonolith)
endif()

if(UTOPIA_INSTALL_LIBMESH)
  include(InstallLIBMESH)
endif()

# ##############################################################################

# #################UTOPIA###################

find_package(Utopia REQUIRED)
if(Utopia_FOUND)
  set(Utopia_FOUND TRUE)
  message(STATUS "Utopia Found.")
  add_definitions(${UTOPIA_DEFS})

  if(NOT UTOPIA_ENABLE_PETSC AND NOT UTOPIA_ENABLE_TRILINOS)
    message(
      FATAL_ERROR
        "Utopia needs to be installed with petsc or trilinos enabled as backends."
    )
  endif()

  list(APPEND UTOPIA_FE_DEP_LIBRARIES ${UTOPIA_LIBRARIES})
  list(APPEND UTOPIA_FE_DEP_INCLUDES ${UTOPIA_INCLUDES})
endif()

# message(STATUS "UTOPIA_FE_DEP_LIBRARIES:${UTOPIA_FE_DEP_LIBRARIES}")

# ##############################################################################
# Temporary fix to target import
# #################YAML-CPP###################
find_package(yaml-cpp HINTS ${UTOPIA_YAML_CPP_DIR} REQUIRED)

# if(UTOPIA_ENABLE_ARBORX) set(ARBORX_SEARCH_PATHS "")

# if(UTOPIA_ENABLE_ENV_READ) set(ARBORX_SEARCH_PATHS
# "${ARBORX_SEARCH_PATHS};$ENV{ARBORX_DIR}") endif()

# find_package(ArborX HINTS ${ARBORX_SEARCH_PATHS} REQUIRED)

# if(ArborX_FOUND) # Add includes to build_includes, .... endif() endif()

if(UTOPIA_ENABLE_SFEM)
  include(FetchContent)

  if(NOT UTOPIA_MATRIX_IO_DIR AND NOT MatrixIO_DIR)
    message(FATAL_ERROR "SFEM required MatrixIO_DIR to be defined!")
  endif()

  set(MatrixIO_DIR ${UTOPIA_MATRIX_IO_DIR}) 
  set(SFEM_ENABLE_SUBMODULE ON)
  set(SFEM_ENABLE_SUBMODULE ON CACHE INTERNAL "")

  FetchContent_Declare(
  sfem
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/sfem
  )

  FetchContent_MakeAvailable(sfem)
  add_library(SFEM::sfem ALIAS sfem)
  list(APPEND UTOPIA_FE_SUBMODULES sfem)
endif()

if(UTOPIA_ENABLE_LIBMESH)
  if(UTOPIA_ENABLE_STK)
    message(
      FATAL_ERROR
        "UtopiaFE cannot be compiled with libmesh and stk enabled at the same time."
    )
  endif()
  if(NOT UTOPIA_INSTALL_LIBMESH)
    find_package(LIBMESH REQUIRED)
  else()
    find_package(LIBMESH QUIET)
  endif()
  if(LIBMESH_FOUND)
    message(STATUS "Libmesh found.")
    list(APPEND UTOPIA_FE_DEP_LIBRARIES ${LIBMESH_LIBRARIES})
    list(APPEND UTOPIA_FE_DEP_INCLUDES ${LIBMESH_INCLUDE_DIR})
  else()
    if(UTOPIA_INSTALL_LIBMESH)
      message(STATUS "[Instructions] Install Libmesh with `export PETSC_DIR=<path_to_petsc>; make libmesh`")
    endif()
  endif()
endif()

# if(UTOPIA_ENABLE_MOONOLITH)
#   if(NOT UTOPIA_INSTALL_MOONOLITH)
#     find_package(ParMoonolith REQUIRED)
#   else()
#     find_package(ParMoonolith QUIET)
#   endif()
#   if(ParMoonolith_FOUND)
#     set(ParMoonolith_FOUND TRUE)
#     message(STATUS "ParMoonolith found.")
#     list(APPEND UTOPIA_FE_DEP_LIBRARIES ${MOONOLITH_LIBRARIES})
#     list(APPEND UTOPIA_FE_DEP_INCLUDES ${MOONOLITH_INCLUDES})
#   endif()
# endif()

if(UTOPIA_ENABLE_MOONOLITH AND NOT UTOPIA_INSTALL_MOONOLITH)
  if(UTOPIA_INSTALL_MOONOLITH_EXTERNAL)
    find_package(ParMoonolith QUIET)
  else()
    find_package(ParMoonolith REQUIRED)
  endif()
    
  set(ParMoonolith_FOUND TRUE)
  message(STATUS "ParMoonolith found.")
  list(APPEND UTOPIA_FE_DEP_LIBRARIES ${MOONOLITH_LIBRARIES})
  list(APPEND UTOPIA_FE_DEP_INCLUDES ${MOONOLITH_INCLUDES})
endif()



if(UTOPIA_ENABLE_INTREPID2)
  set(UTOPIA_INTREPID2_SEARCH_PATHS "${UTOPIA_TRILINOS_DIR}/../Intrepid2")
  find_package(Intrepid2 HINTS ${UTOPIA_INTREPID2_SEARCH_PATHS} REQUIRED)
  if(Intrepid2_FOUND)
    message(STATUS "Intrepid2 found.")
    if(Trilinos_Kokkos_FOUND)
      list(APPEND UTOPIA_FE_DEP_LIBRARIES ${Intrepid2_LIBRARIES})
      list(APPEND UTOPIA_FE_DEP_INCLUDES ${Trilinos_INCLUDE_DIRS})
      # message(STATUS "Trilinos_INCLUDE_DIRS;${Trilinos_INCLUDE_DIRS}")
    else()
      message(
        FATAL_ERROR "UtopiaFE needs kokkos installed with utopia trilinos.")
    endif()
  endif()
endif()

if(UTOPIA_ENABLE_STK)
  if(UTOPIA_ENABLE_LIBMESH)
    message(
      FATAL_ERROR
        "UtopiaFE cannot be compiled with libmesh and stk enabled at the same time."
    )
  endif()
  set(UTOPIA_STK_SEARCH_PATHS "${Trilinos_DIR}/../STK")
  if(UTOPIA_ENABLE_ENV_READ)
    set(UTOPIA_STK_SEARCH_PATHS
        "${UTOPIA_STK_SEARCH_PATHS};$ENV{TRILINOS_DIR}/lib/cmake/STK")
  endif()
  find_package(STK HINTS ${UTOPIA_STK_SEARCH_PATHS} REQUIRED)
  if(STK_FOUND)
	  message(STATUS "STK found.")
    list(APPEND UTOPIA_FE_DEP_LIBRARIES ${STK_LIBRARIES})
    list(APPEND UTOPIA_FE_DEP_INCLUDES ${Trilinos_INCLUDE_DIRS})
  endif()
endif()

if(UTOPIA_ENABLE_MARS)
  find_package(Mars REQUIRED)
  if(Mars_FOUND)
    set(Mars_FOUND TRUE)
    message(STATUS "Mars found!")
    get_target_property(MARS_DEP_LIBRARIES Mars::mars INTERFACE_LINK_LIBRARIES)
    list(APPEND UTOPIA_FE_DEP_LIBRARIES ${MARS_LIBRARIES} ${MARS_DEP_LIBRARIES})
    list(APPEND UTOPIA_FE_DEP_INCLUDES ${MARS_INCLUDES})
  else()
    message(FATAL_ERROR "MARS not found!")
  endif()
endif()

# ##############################################################################

# message(STATUS "UTOPIA_DIR: ${UTOPIA_DIR}")


macro(print_dependency_table)

  set(SMALL_DEP_TABLE
      "\n_______________________________\n\n   BACKENDS and STATUS TABLE\n")
  set(SMALL_DEP_TABLE "${SMALL_DEP_TABLE}-------------------------------\n")
  set(SMALL_DEP_TABLE
      "${SMALL_DEP_TABLE}backend\t|status\t|found\n-------------------------------\n"
  )
  set(SMALL_DEP_TABLE
      "-${SMALL_DEP_TABLE}mpi\t|${UTOPIA_ENABLE_MPI}\t|${MPI_FOUND}\n")

  set(SMALL_DEP_TABLE
      "-${SMALL_DEP_TABLE}Utopia\t|ON\t|${Utopia_FOUND}\n")

  set(SMALL_DEP_TABLE
      "-${SMALL_DEP_TABLE}Libmesh\t|${UTOPIA_ENABLE_LIBMESH}\t|${LIBMESH_FOUND}\n"
  )

  set(SMALL_DEP_TABLE
      "-${SMALL_DEP_TABLE}ParMoon\t|${UTOPIA_ENABLE_MOONOLITH}\t|${ParMoonolith_FOUND}\n"
  )

  set(SMALL_DEP_TABLE
      "-${SMALL_DEP_TABLE}Intrep\t|${UTOPIA_ENABLE_INTREPID2}\t|${INTREPID_2_FOUND}\n"
  )

  set(SMALL_DEP_TABLE
      "-${SMALL_DEP_TABLE}Mars\t|${UTOPIA_ENABLE_MARS}\t|${Mars_FOUND}\n")

  set(SMALL_DEP_TABLE
      "-${SMALL_DEP_TABLE}SFEM\t|${UTOPIA_ENABLE_SFEM}\t| --\n")

  set(SMALL_DEP_TABLE "${SMALL_DEP_TABLE}_______________________________\n")

  message(STATUS ${SMALL_DEP_TABLE})
endmacro()

macro(log_dependency_table)

  set(DEP_TABLE "backends:\n")
  set(DEP_TABLE "${DEP_TABLE}  - mpi:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_MPI: ${UTOPIA_ENABLE_MPI}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MPI_DIR: ${UTOPIA_MPI_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MPI_VERSION: ${UTOPIA_MPI_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - Utopia:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_DIR: ${UTOPIA_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_VERSION: ${UTOPIA_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - Libmesh:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_LIBMESH: ${UTOPIA_ENABLE_LIBMESH}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_LIBMESH_DIR: ${UTOPIA_LIBMESH_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_LIBMESH_VERSION: ${UTOPIA_LIBMESH_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - Moonolith:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_MOONOLITH: ${UTOPIA_ENABLE_MOONOLITH}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MOONOLITH_DIR: ${UTOPIA_MOONOLITH_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MOONOLITH_VERSION: ${UTOPIA_MOONOLITH_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - Intrepid2:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_INTREPID2: ${UTOPIA_ENABLE_INTREPID2}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_INTREPID2_DIR: ${UTOPIA_INTREPID2_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_INTREPID2_VERSION: ${UTOPIA_INTREPID2_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - STK:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_STK: ${UTOPIA_ENABLE_STK}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_STK_DIR: ${UTOPIA_STK_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_STK_VERSION: ${UTOPIA_STK_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - Mars:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_MARS: ${UTOPIA_ENABLE_MARS}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MARS_DIR: ${UTOPIA_MARS_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MARS_VERSION: ${UTOPIA_MARS_VERSION}\n")

  # OPTIONS

  set(DEP_TABLE "${DEP_TABLE}options:\n")
  getlistofvarsstartingwith("UTOPIA_" matchedVars)
  foreach(_var IN LISTS matchedVars)
    set(DEP_TABLE "${DEP_TABLE}  - ${_var}: ${${_var}}\n")
  endforeach()

  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfig.yaml" ${DEP_TABLE})
endmacro()
