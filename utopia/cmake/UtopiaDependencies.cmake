# Cmake for checking which dependencies are enabled and if so then point to
# correct cmake subdirectory.

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE)
  message(
    FATAL_ERROR
      "In-source builds are not allowed.
        Please create a directory and run cmake from there, passing the path
        to this source directory as the last argument.
        This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
        Please delete them.")
endif()

if(CYGWIN)
  include(InstallPetscCygwin)
endif()

if(LINUX)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

if(USE_SPIKE_SOLVERS)
  list(APPEND UTOPIA_MODULES spike/solvers)
endif()

if(UTOPIA_ENABLE_BLAS)
  add_subdirectory(backend/blas)
endif()

# if(UTOPIA_ENABLE_PASSO_EXTENSIONS) set(WITH_PASSO_EXTENSIONS TRUE) endif()

if(UTOPIA_ENABLE_CXX14_FEATURES)
  set(UTOPIA_WITH_CPP14 TRUE)
endif()

if(UTOPIA_ENABLE_DEPRECATED_API)
  set(UTOPIA_DEPRECATED_API ON)
endif()

if(UTOPIA_ENABLE_GPERFTOOLS)
  find_package(Gperftools)

  if(Gperftools_FOUND)
    link_libraries(gperftools::profiler)
  else()
    message(WARNING "GPERFTOOLS NOT FOUND")
  endif()
endif()

if(UTOPIA_ENABLE_LOCK_CHECKING)
  set(UTOPIA_ENABLE_LOCK_CHECK TRUE)
endif()

if(UTOPIA_ENABLE_METIS)
  add_subdirectory(backend/metis)
endif()

if(UTOPIA_ENABLE_PARMETIS)
  add_subdirectory(backend/parmetis)
endif()

# using this to defining UTOPIA_LAMBDA correctly
if(UTOPIA_ENABLE_PETSC)
  add_subdirectory(backend/petsc)
endif()

if(UTOPIA_ENABLE_SCRIPTING)
  add_subdirectory(scripting)
endif()

# backend modules
if(UTOPIA_ENABLE_TRACE)
  set(UTOPIA_TRACE_ENABLED ON)
endif()

if(UTOPIA_ENABLE_TRACE_EXPR)
  set(UTOPIA_ENABLE_TRACE ON)
  set(UTOPIA_TRACE_EXPR_ENABLED ON)
  set(UTOPIA_TRACE_ENABLED ON)
endif()

if(UTOPIA_ENABLE_TRILINOS)
  add_subdirectory(backend/trilinos)
endif()

if(UTOPIA_ENABLE_TRILINOS)
  find_package(Trilinos)
  if(Trilinos_FOUND)
    # include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS}
    # ${Trilinos_TPL_INCLUDE_DIRS}) set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}
    # ${Kokkos_CXX_FLAGS}") message(STATUS
    # "Kokkos_CXX_FLAGS=${Kokkos_CXX_FLAGS},${CMAKE_CXX_FLAGS}") message(STATUS
    # "HERE:
    # ${Trilinos_CXX_COMPILER_FLAGS},${Kokkos_CXX_FLAGS},${TRILINOS_DEFINITIONS},${CMAKE_CXX_STANDARD},${INTERFACE_COMPILE_FEATURES}")
    set(UTOPIA_WITH_TRILINOS ON)
  endif()
endif()

if(UTOPIA_ENABLE_VC)
  add_subdirectory(backend/vc)
endif()

if(UTOPIA_ENABLE_YAML_CPP)
  add_subdirectory(backend/yamlcpp)
endif()

if(UTOPIA_INSTALL_PETSC AND UTOPIA_ENABLE_CYGWIN)
  include(InstallPetsc)

else()
  include(InstallPetsc)
endif()

if(UTOPIA_INSTALL_PETSC_DEBUG)
  include(InstallPetscDebug)
endif()

if(UTOPIA_ENABLE_POLYMORPHIC)
  add_subdirectory(backend/polymorphic)
endif()
if(UTOPIA_INSTALL_TRILINOS)
  include(InstallTrilinos)
endif()

# weird stuff goes here
if(UTOPIA_STATIC_DEPENDENCIES_ONLY)
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  endif()
endif()

message(STATUS "[Status] UTOPIA_ROOT_PATH: ${UTOPIA_ROOT_PATH}")
include(${UTOPIA_ROOT_PATH}/cmake/UtopiaCompilerFeatures.cmake)

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################MPI######################
find_package(MPIExtended)
if(UTOPIA_ENABLE_MPI)
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
    endif()

    set(UTOPIA_WITH_MPI ON)
    set(UTOPIA_MPI TRUE)

    if(!MPI_DIR)
      set(MPI_DIR ${MPI_DIR})
    endif()
  endif()
else()
  message(WARNING "NO Proper MPI installation")
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################DOXYGEN######################

find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  message(
    WARNING
      "Doxygen is needed to build the documentation. Please install it correctly"
  )
endif()
# -- Configure the Template Doxyfile for our specific project
configure_file(${UTOPIA_ROOT_PATH}/Doxyfile.txt ${CMAKE_BINARY_DIR} @ONLY
               IMMEDIATE)
# -- Add a custom target to run Doxygen when ever the project is built
add_custom_target(
  docs
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile.txt
  SOURCES ${CMAKE_BINARY_DIR}/Doxyfile.txt)
# IF you do NOT want the documentation to be generated EVERY time you build the
# project then leave out the 'ALL' keyword from the above command.
