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

# #################  BACKENDS  ######################

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

    set(UTOPIA_ENABLE_MPI ON)
    set(UTOPIA_ENABLE_MPI TRUE)

    set(UTOPIA_MPI_DIR ${MPI_LIBRARIES})
    set(UTOPIA_MPI_VERSION ${MPI_C_VERSION})
  endif()
else()
  message(WARNING "NO Proper MPI installation")
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################BLAS, LAPACK, UMFPACK######################

if(UTOPIA_ENABLE_BLAS)
  # find dependencies ###################################################
  if(PITZ_DORA)
    find_package(OpenBLAS REQUIRED)
    if(OPEN_BLAS_FOUND)
      list(APPEND UTOPIA_BUILD_INCLUDES ${BLAS_INCLUDE_DIR})
      list(APPEND UTOPIA_DEP_LIBRARIES ${BLAS_LIBRARIES})
      list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
      set(UTOPIA_ENABLE_BLAS TRUE)
      # set(UTOPIA_ENABLE_OPEN_BLAS ON) set(UTOPIA_BLAS_DIR ${BLAS_blas_LIBRARY}
      # PARENT_SCOPE) # BLAS_blas_LIBRARY set(UTOPIA_BLAS_DIR ${BLAS_VERSION}
      # PARENT_SCOPE)
    else()
      message(WARNING "[Warning] open-blas not found")
    endif()
  else()
    set(OPEN_BLAS_FOUND FALSE)
    find_package(BLAS REQUIRED)
    if(BLAS_FOUND)
      list(APPEND UTOPIA_BUILD_INCLUDES ${BLAS_INCLUDE_DIR})
      list(APPEND UTOPIA_DEP_LIBRARIES ${BLAS_LIBRARIES})
      list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
      set(UTOPIA_ENABLE_BLAS TRUE)
    else()
      message(WARNING "[Warning] blas not found")
      # set(BLAS_FOUND FALSE)
    endif()
  endif()

  find_package(LAPACK)
  if(LAPACK_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${LAPACK_INCLUDE_DIR})
    list(APPEND UTOPIA_DEP_LIBRARIES ${LAPACK_LIBRARIES})
    list(APPEND UTOPIA_DEFS ${LAPACK_DEFINITIONS})
    # set(UTOPIA_ENABLE_LAPACK ON) set(UTOPIA_ENABLE_LAPACK TRUE)
  else()
    message(WARNING "[Warning] lapack not found")
    # set(UTOPIA_ENABLE_LAPACK OFF)
  endif()

  find_package(UMFPACK)
  if(UMFPACK_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${UMFPACK_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${UMFPACK_LIBRARIES})
    # set(UTOPIA_ENABLE_UMFPACK ON) set(UTOPIA_ENABLE_UMFPACK TRUE PARENT_SCOPE)
  else()
    message(WARNING "[Warning] Umfpack not found")
    # set(UTOPIA_ENABLE_UMFPACK OFF)
  endif()

  add_subdirectory(backend/blas)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################METIS######################

if(UTOPIA_ENABLE_METIS)
  find_package(METIS REQUIRED)
  if(METIS_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${METIS_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${METIS_LIBRARIES})

  else()
    message(WARNING "[Warning] Metis not found")
  endif()
  add_subdirectory(backend/metis)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################PARMETIS######################

if(UTOPIA_ENABLE_PARMETIS)
  find_package(PARMETIS REQUIRED)

  if(NOT PARMETIS_FOUND)
    message(FATAL_ERROR "Unable to find ParMetis")
    # include(FetchContent) message(STATUS "Fetching parmetis, since it could
    # not be found.")

    # if(UTOPIA_DEPENDENCIES_DIR) set(PARMETIS_INSTALL_DIR
    # ${UTOPIA_DEPENDENCIES_DIR}/PARMETIS) else() set(PARMETIS_INSTALL_DIR
    # ${CMAKE_SOURCE_DIR}/../external/PARMETIS) endif()

    # set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage") set(PARMETIS_URL
    # "https://github.com/KarypisLab/ParMETIS")

    # FetchContent_Declare( parmetis GIT_REPOSITORY ${PARMETIS_URL} # GIT_TAG
    # origin/sampler GIT_TAG origin/development )

    # FetchContent_MakeAvailable(parmetis)

    # add_library(ParMetis::parmetis ALIAS parmetis)
    # target_link_libraries(utopia PUBLIC parmetis)

  else()
    # target_include_directories(utopia PUBLIC ${PARMETIS_INCLUDES})
    # target_link_libraries(utopia PUBLIC ${PARMETIS_LIBRARIES})

    list(APPEND UTOPIA_BUILD_INCLUDES ${PARMETIS_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${PARMETIS_LIBRARIES})

  endif()
  add_subdirectory(backend/parmetis)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################PETSC######################
if(UTOPIA_ENABLE_PETSC)

  set(PETSC_TEST_RUNS TRUE)
  set(PETSC_EXECUTABLE_RUNS TRUE) # On daint we cannot run them

  find_package(PETSc REQUIRED)

  if(PETSC_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${PETSC_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${PETSC_LIBRARIES})

    set(CMAKE_C_COMPILER ${PETSC_COMPILER})

    if(NOT MPI_CXX_COMPILER)
      set(MPI_CXX_COMPILER $ENV{MPI_CXX_COMPILER})
      message(STATUS "compiler ${MPI_CXX_COMPILER}")
    endif()

    if(MPI_CXX_COMPILER)
      set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
      set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
    else()
      execute_process(COMMAND mpicxx -v RESULT_VARIABLE MPICXX_FAILED)

      if(MPICXX_FAILED)
        message(
          STATUS
            "Using CMAKE compiler, you can define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>"
        )
      else()
        message(STATUS "-----------------------------------------------")
        message(
          STATUS
            "\n[MPI] using mpicxx for compiling c++ files.\nIf you want to use your own compiler define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>"
        )
        message(STATUS "-----------------------------------------------")
        set(CMAKE_CXX_COMPILER mpicxx)
        set(CMAKE_CXX_COMPILER_DEBUG mpicxx)
      endif()
    endif()

    # set(UTOPIA_ENABLE_PETSC ON)

    set(UTOPIA_PETSC_DIR ${PETSC_DIR})
    set(UTOPIA_PETSC_VERSION ${PETSC_VERSION})

  else()
    message(WARNING "[Warning] Petsc not found")
    set(UTOPIA_ENABLE_PETSC FALSE)
  endif()

  if(PETSC_FOUND AND UTOPIA_ENABLE_SLEPC)
    find_package(SLEPC REQUIRED)
    if(SLEPC_FOUND)
      list(APPEND UTOPIA_BUILD_INCLUDES ${SLEPC_INCLUDES})
      list(APPEND UTOPIA_DEP_LIBRARIES ${SLEPC_LIBRARIES})
      message(STATUS "Slepc FOUND")

      set(UTOPIA_SLEPC_DIR ${SLEPC_DIR})
      set(UTOPIA_SLEPC_VERSION ${SLEPC_VERSION})
      # set(UTOPIA_ENABLE_SLEPC TRUE PARENT_SCOPE)
    else()
      message(WARNING "[Warning] Slepc not found")
      # set(UTOPIA_ENABLE_SLEPC FALSE)
    endif()
  endif()
  add_subdirectory(backend/petsc)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################TRILINOS######################

if(UTOPIA_ENABLE_TRILINOS)
  # find dependencies
  find_package(Trilinos REQUIRED)
  if(Trilinos_FOUND)
    # These lines are needed so the utopia-config.makefile will be correctly
    # built. The UtopiaTargets.cmake file is populated automatically by cmake
    # ##########################################################################
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
      # e.g Amesos2: UTOPIA_WITH_TRILINOS_AMESOS2=TRUE if package is found in
      # Trilinos_PACKAGE_LIST
      string(TOUPPER ${package} packageUpper)
      string(TOLOWER ${package} packageLower)
      list(FIND Trilinos_PACKAGE_LIST ${package} PACKAGE_FOUND)
      if(NOT PACKAGE_FOUND EQUAL -1)
        set(UTOPIA_WITH_TRILINOS_${packageUpper} TRUE)
        if(${package} STREQUAL "MueLu")
          # ugly hack, but we need to link with muelu-adapters also I cannot
          # wait until Trilinos finally supports cmake targets correctly
          list(APPEND UTOPIA_TRILINOS_DEPS "muelu-adapters")
        endif()
        list(APPEND UTOPIA_TRILINOS_DEPS ${packageLower})
      endif()
    endforeach()

    find_package(TpetraExt)
    if(TRILINOS_TPETRAEXT_FOUND)
      set(UTOPIA_WITH_TRILINOS_TPETRAEXT TRUE)
    endif()
  else()
    message(WARNING "[Warning] Trilinos not found")
  endif()
  add_subdirectory(backend/trilinos)
endif()

# if(UTOPIA_ENABLE_TRILINOS)

# find_package(Trilinos) # if(Trilinos_FOUND) #   # include_directories(SYSTEM
# ${Trilinos_INCLUDE_DIRS} #   # ${Trilinos_TPL_INCLUDE_DIRS}) #   #
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} # # ${Kokkos_CXX_FLAGS}") #   #
# message(STATUS "Kokkos_CXX_FLAGS=${Kokkos_CXX_FLAGS},${CMAKE_CXX_FLAGS}") # #
# message( # #   STATUS #   #     "HERE: #   #
# ${Trilinos_CXX_COMPILER_FLAGS},${Kokkos_CXX_FLAGS},${TRILINOS_DEFINITIONS},${CMAKE_CXX_STANDARD},${INTERFACE_COMPILE_FEATURES}"
# #   # ) #   # set(UTOPIA_ENABLE_TRILINOS ON) # endif() endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################YAML######################

if(UTOPIA_ENABLE_YAML_CPP)
  find_package(
    yaml-cpp HINTS ${CMAKE_SOURCE_DIR}/../../dependencies/
    ${CMAKE_SOURCE_DIR}/../../dependencies/yaml-cpp/share/cmake/yaml-cpp)

  if(yaml-cpp_FOUND)
    set(UTOPIA_ENABLE_YAML_CPP ON)
    set(UTOPIA_ENABLE_YAML_CPP TRUE)

    set(UTOPIA_YAML_CPP_DIR ${yaml-cpp_DIR})

    list(APPEND YAMLCPP_MODULES .)

    # utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${YAMLCPP_MODULES}")

    get_target_property(YAML_CPP_INCLUDE_DIR yaml-cpp
                        INTERFACE_INCLUDE_DIRECTORIES)

    get_target_property(YAML_CPP_LIBRARIES yaml-cpp IMPORTED_LOCATION)

    if(NOT YAML_CPP_LIBRARIES)
      get_target_property(YAML_CPP_LIBRARIES yaml-cpp IMPORTED_LOCATION_RELEASE)
    endif()
    if(NOT YAML_CPP_LIBRARIES)
      get_target_property(YAML_CPP_LIBRARIES yaml-cpp IMPORTED_LOCATION_DEBUG)
    endif()
    if(NOT YAML_CPP_LIBRARIES)
      get_target_property(YAML_CPP_LIBRARIES yaml-cpp
                          IMPORTED_LOCATION_NOCONFIG)
    endif()

    message(
      STATUS
        "yaml-cpp found! Includes: ${YAML_CPP_INCLUDE_DIR}\nLibrary: ${YAML_CPP_LIBRARIES}\n"
    )

    # target_link_libraries(utopia PUBLIC yaml-cpp)
    # target_include_directories(utopia PUBLIC ${YAML_CPP_INCLUDE_DIR})
    # target_link_libraries(utopia PUBLIC ${YAML_CPP_LIBRARIES})

    list(APPEND UTOPIA_BUILD_INCLUDES ${YAML_CPP_INCLUDE_DIR})
    list(APPEND UTOPIA_DEP_LIBRARIES ${YAML_CPP_LIBRARIES})

    set(UTOPIA_BUILD_INCLUDES ${UTOPIA_BUILD_INCLUDES})

    set(UTOPIA_DEP_LIBRARIES ${UTOPIA_DEP_LIBRARIES})

  else()
    include(../../cmake/InstallYAMLCPP.cmake)
  endif()
  add_subdirectory(backend/yamlcpp)
endif()

if(UTOPIA_ENABLE_VC)
  add_subdirectory(backend/vc)
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

# if(UTOPIA_ENABLE_PASSO_EXTENSIONS) set(WITH_PASSO_EXTENSIONS TRUE) endif()

if(UTOPIA_ENABLE_CXX14_FEATURES)
  set(UTOPIA_ENABLE_CPP14 TRUE)
endif()

if(UTOPIA_ENABLE_CXX17_FEATURES)
  set(UTOPIA_ENABLE_CPP17 TRUE)
  set(UTOPIA_ENABLE_CPP14 TRUE)
endif()

# If both enabled just use if(UTOPIA_ENABLE_CXX14_FEATURES AND
# UTOPIA_ENABLE_CXX17_FEATURES) set(UTOPIA_ENABLE_CPP14 FALSE)
# set(UTOPIA_ENABLE_CPP17 TRUE) endif()

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

# weird stuff goes here
if(UTOPIA_STATIC_DEPENDENCIES_ONLY)
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  endif()
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

message(STATUS "[Status] UTOPIA_ROOT_PATH: ${UTOPIA_ROOT_PATH}")
include(${UTOPIA_ROOT_PATH}/cmake/UtopiaCompilerFeatures.cmake)

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
