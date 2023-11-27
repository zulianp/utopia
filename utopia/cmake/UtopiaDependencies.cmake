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

if(UTOPIA_ENABLE_VC)
  add_subdirectory(backend/vc)
endif()

if(UTOPIA_ENABLE_POLYMORPHIC)
  add_subdirectory(backend/polymorphic)
endif()

if(UTOPIA_ENABLE_CXX14_FEATURES)
  set(UTOPIA_ENABLE_CPP14 TRUE)
endif()

if(UTOPIA_ENABLE_CXX17_FEATURES)
  set(UTOPIA_ENABLE_CPP17 TRUE)
  set(UTOPIA_ENABLE_CPP14 TRUE)
endif()

if(UTOPIA_ENABLE_DEPRECATED_API)
  set(UTOPIA_DEPRECATED_API ON)
endif()

# Static dependencies.
if(UTOPIA_STATIC_DEPENDENCIES_ONLY)
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  endif()
endif()

# #################  BACKENDS  ######################

# ##############################################################################
# ##############################################################################
# ##############################################################################

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

# ##############################################################################
# ##############################################################################
# ##############################################################################

# ####BLAS, LAPACK, UMFPACK####

if(UTOPIA_ENABLE_BLAS)
  find_package(BLAS REQUIRED)
  if(BLAS_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${BLAS_INCLUDE_DIR})
    list(APPEND UTOPIA_DEP_LIBRARIES ${BLAS_LIBRARIES})
    list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})

    set(UTOPIA_ENABLE_BLAS ON)
    set(UTOPIA_BLAS_DIR ${BLAS_blas_LIBRARY})

  else()
    message(WARNING "[Warning] Blas not found")
  endif()

  find_package(LAPACK)
  if(LAPACK_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${LAPACK_INCLUDE_DIR})
    list(APPEND UTOPIA_DEP_LIBRARIES ${LAPACK_LIBRARIES})
    list(APPEND UTOPIA_DEFS ${LAPACK_DEFINITIONS})

    set(UTOPIA_LAPACK_DIR ${LAPACK_LIBRARIES})
  else()
    message(WARNING "[Warning] Lapack not found")
  endif()

  find_package(UMFPACK)
  if(UMFPACK_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${UMFPACK_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${UMFPACK_LIBRARIES})

    set(UTOPIA_UMFPACK_DIR ${UMFPACK_LIBDIR})
  else()
    message(WARNING "[Warning] Umfpack not found")
  endif()

  add_subdirectory(backend/blas)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################METIS####################

if(UTOPIA_ENABLE_METIS)
  find_package(METIS REQUIRED)
  if(METIS_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${METIS_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${METIS_LIBRARIES})
    set(UTOPIA_ENABLE_METIS ON)
  else()
    message(WARNING "[Warning] Metis not found")
  endif()
  add_subdirectory(backend/metis)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################PARMETIS#################

if(UTOPIA_ENABLE_PARMETIS)
  find_package(ParMetis REQUIRED)

  if(ParMetis_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${ParMetis_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${ParMetis_LIBRARIES})
  else()
    message(FATAL_ERROR "[Warning] ParMetis not found.")
  endif()
  add_subdirectory(backend/parmetis)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################
# ####LOCAL_INSTALL_OPTIONS####

if(UTOPIA_INSTALL_PETSC
   AND NOT CYGWIN
   AND UTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL)
  include(InstallPetsc)
endif()

if(UTOPIA_INSTALL_TRILINOS AND UTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL)
  include(InstallTrilinos)
endif()

if(UTOPIA_INSTALL_SLEPC AND UTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL)
  include(InstallSlepc)
endif()

if(UTOPIA_INSTALL_YAML_CPP)
  include(InstallYAMLCPP)
endif() 

# #################PETSC####################
if(UTOPIA_ENABLE_PETSC)

  set(PETSC_TEST_RUNS TRUE)
  set(PETSC_EXECUTABLE_RUNS TRUE) # On daint we cannot run them

  if(NOT UTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL)
    find_package(PETSc REQUIRED)
  else()
    find_package(PETSc QUIET)
  endif()

  if(PETSC_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${PETSC_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${PETSC_LIBRARIES})
    list(APPEND UTOPIA_DEP_INCLUDES ${PETSC_INCLUDES})

    set(CMAKE_C_COMPILER ${PETSC_COMPILER})

    if(NOT MPI_CXX_COMPILER)
      set(MPI_CXX_COMPILER $ENV{MPI_CXX_COMPILER})
      message(STATUS "compiler ${MPI_CXX_COMPILER}")
    endif()

    if(NOT MPI_CXX_COMPILER)
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

    set(UTOPIA_PETSC_DIR ${PETSC_DIR})
    set(UTOPIA_PETSC_VERSION ${PETSC_VERSION})

  else()
    message(WARNING "[Warning] Petsc not found")
  endif()

  if(PETSC_FOUND AND UTOPIA_ENABLE_SLEPC)
    find_package(SLEPc QUIET)
    if(SLEPC_FOUND)
      list(APPEND UTOPIA_BUILD_INCLUDES ${SLEPC_INCLUDES})
      list(APPEND UTOPIA_DEP_LIBRARIES ${SLEPC_LIBRARIES})
      message(STATUS "Slepc FOUND")

      set(UTOPIA_SLEPC_DIR ${SLEPC_DIR})
      set(UTOPIA_SLEPC_VERSION ${PETSC_VERSION})
    else()
      message(WARNING "[Warning] Slepc not found")
    endif()
  endif()
  add_subdirectory(backend/petsc)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################TRILINOS#################

if(UTOPIA_ENABLE_TRILINOS)

  list(APPEND Trilinos_SEARCH_PATHS "${Trilinos_DIR};${TRILINOS_DIR}")

  if(UTOPIA_ENABLE_ENV_READ)
    list(APPEND Trilinos_SEARCH_PATHS "$ENV{Trilinos_DIR};$ENV{TRILINOS_DIR}")
  endif()

  set(Trilinos_SEARCH_PATHS
      "${Trilinos_SEARCH_PATHS};${TRILINOS_DIR}/lib/cmake/Trilinos;${Trilinos_DIR}/lib/cmake/Trilinos;${TRILINOS_DIR}/lib64/cmake/Trilinos;${Trilinos_DIR}/lib64/cmake/Trilinos"
  )

  # find dependencies
  if(NOT UTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL)
    set(Trilinos_FOUND FALSE)
    find_package(Trilinos PATHS ${Trilinos_SEARCH_PATHS} REQUIRED)
  else()
    set(Trilinos_FOUND FALSE)
    set(Trilinos_SEARCH_PATHS
        "${CMAKE_SOURCE_DIR}/../external/Trilinos/lib/cmake/Trilinos")
    if(NOT APPLE AND NOT WIN32)
      set(Trilinos_SEARCH_PATHS
          "${CMAKE_SOURCE_DIR}/../external/Trilinos/lib64/cmake/Trilinos")
    endif()
    find_package(Trilinos PATHS ${Trilinos_SEARCH_PATHS} NO_DEFAULT_PATH)
  endif()
  if(Trilinos_FOUND)
    set(Trilinos_FOUND TRUE)
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

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################YAML######################

if(UTOPIA_ENABLE_YAML_CPP)

  set(YAML_CPP_SEARCH_PATHS "${YAMLCPP_INSTALL_DIR};${CMAKE_SOURCE_DIR}/../external/yaml-cpp/")
  if(UTOPIA_ENABLE_ENV_READ)
    set(YAML_CPP_SEARCH_PATHS "${YAML_CPP_SEARCH_PATHS};$ENV{YAMLCPP_DIR}")
  endif()

  find_package(yaml-cpp HINTS ${YAML_CPP_SEARCH_PATHS})

  if(yaml-cpp_FOUND)
    set(yaml-cpp_FOUND TRUE)
    set(UTOPIA_ENABLE_YAML_CPP ON)

    set(UTOPIA_YAML_CPP_DIR ${yaml-cpp_DIR})
    set(UTOPIA_YAML_CPP_VERSION ${yaml-cpp_VERSION})

    list(APPEND YAMLCPP_MODULES .)

    if(APPLE)

      get_target_property(YAML_CPP_INCLUDE_DIR yaml-cpp::yaml-cpp
                          INTERFACE_INCLUDE_DIRECTORIES)

      get_target_property(YAML_CPP_LIBRARIES yaml-cpp::yaml-cpp
                          IMPORTED_LOCATION)

      if(NOT YAML_CPP_LIBRARIES)
        get_target_property(YAML_CPP_LIBRARIES yaml-cpp::yaml-cpp
                            IMPORTED_LOCATION_RELEASE)
      endif()
      if(NOT YAML_CPP_LIBRARIES)
        get_target_property(YAML_CPP_LIBRARIES yaml-cpp::yaml-cpp
                            IMPORTED_LOCATION_DEBUG)
      endif()
      if(NOT YAML_CPP_LIBRARIES)
        get_target_property(YAML_CPP_LIBRARIES yaml-cpp::yaml-cpp
                            IMPORTED_LOCATION_NOCONFIG)
      endif()
    endif()

    message(STATUS "yaml-cpp found!")
    message(STATUS "Includes: ${YAML_CPP_INCLUDE_DIR}")
    message(STATUS "Library: ${YAML_CPP_LIBRARIES}")
    
    list(APPEND UTOPIA_BUILD_INCLUDES ${YAML_CPP_INCLUDE_DIR})
    list(APPEND UTOPIA_DEP_LIBRARIES ${YAML_CPP_LIBRARIES})

    set(UTOPIA_BUILD_INCLUDES ${UTOPIA_BUILD_INCLUDES})
    set(UTOPIA_DEP_LIBRARIES ${UTOPIA_DEP_LIBRARIES})
    add_subdirectory(backend/yamlcpp)
  else()
    message(
      WARNING
        "Help message:\n"
        "---------------------------------------------------------------\n"
        "yaml-cpp not found! To install locally in UTOPIA_DEPENDENCIES_DIR then run `make yaml-cpp` and re-run cmake with options `-DYAMLCPP_DIR=${YAMLCPP_INSTALL_DIR}/lib/cmake/yaml-cpp. Otherwise export YAMLCPP_DIR.`\n"
        "---------------------------------------------------------------\n")
  endif()
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################GPERF-TOOLS##############

if(UTOPIA_ENABLE_GPERFTOOLS)
  find_package(Gperftools)

  if(Gperftools_FOUND)
    link_libraries(gperftools::profiler)
  else()
    message(WARNING "GPERFTOOLS NOT FOUND")
  endif()
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################DOXYGEN##################

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

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################SWIG#####################
if(UTOPIA_ENABLE_SCRIPTING)
  find_package(SWIG REQUIRED)

  if(SWIG_FOUND)
    include(${SWIG_USE_FILE})

    find_package(
      Python3
      COMPONENTS Development NumPy
      REQUIRED)

    if(Python3_FOUND)
      set_property(SOURCE scripting/utopia.i PROPERTY CPLUSPLUS ON)
      set_source_files_properties(scripting/utopia.i PROPERTIES SWIG_FLAGS
                                                                "-includeall")

      set_source_files_properties(SOURCE scripting/utopia.i PROPERTY
                                  SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)

      message(STATUS "Found python libraries: ${Python_LIBRARIES}")
      swig_add_library(
        utopya
        TYPE SHARED
        LANGUAGE python
        SOURCES scripting/utopia.i scripting/utopia_script.hpp
                scripting/utopia_script.cpp)

      # swig_link_libraries(utopya PUBLIC Python3::Python Python3::NumPy utopia)
      swig_link_libraries(utopya PUBLIC Python3::Module utopia)
      set_target_properties(utopya PROPERTIES SUFFIX ".so")

      set_property(TARGET utopya PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES
                                          TRUE)

      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia.py DESTINATION lib)

      install(TARGETS ${SWIG_MODULE_utopya_REAL_NAME} LIBRARY DESTINATION lib)

      # find_package(NumPy REQUIRED)
      message(STATUS "PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")

    else()
      message(WARNING "Could not find Python - Disabling scripting support")
    endif()
    add_subdirectory(scripting/)
  else()
    message(WARNING "Could not find SWIG - Disabling scripting support")
  endif()
endif()

# ##############################################################################

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
      "${SMALL_DEP_TABLE}petsc\t|${UTOPIA_ENABLE_PETSC}\t|${PETSC_FOUND}\n")
  set(SMALL_DEP_TABLE
      "${SMALL_DEP_TABLE}slepc\t|${UTOPIA_ENABLE_SLEPC}\t|${SLEPC_FOUND}\n")
  set(SMALL_DEP_TABLE
      "${SMALL_DEP_TABLE}trilinos|${UTOPIA_ENABLE_TRILINOS}\t|${Trilinos_FOUND}\n"
  )
  set(SMALL_DEP_TABLE
      "${SMALL_DEP_TABLE}blas\t|${UTOPIA_ENABLE_BLAS}\t|${BLAS_FOUND}\n")
  set(SMALL_DEP_TABLE
      "${SMALL_DEP_TABLE}lapack\t|${UTOPIA_ENABLE_PETSC}\t|${LAPACK_FOUND}\n")
  set(SMALL_DEP_TABLE
      "${SMALL_DEP_TABLE}umfpack\t|${UTOPIA_ENABLE_PETSC}\t|${UMFPACK_FOUND}\n")
  set(SMALL_DEP_TABLE
      "${SMALL_DEP_TABLE}yaml\t|${UTOPIA_ENABLE_YAML_CPP}\t|${yaml-cpp_FOUND}\n"
  )
  set(SMALL_DEP_TABLE "${SMALL_DEP_TABLE}_______________________________\n")

  message(STATUS ${SMALL_DEP_TABLE})
endmacro()

macro(log_dependency_table)

  set(DEP_TABLE "backends:\n")
  set(DEP_TABLE "${DEP_TABLE}  - mpi:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_MPI: ${UTOPIA_ENABLE_MPI}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MPI_DIR: ${UTOPIA_MPI_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_MPI_VERSION: ${UTOPIA_MPI_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - petsc:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_PETSC: ${UTOPIA_ENABLE_PETSC}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_PETSC_DIR: ${UTOPIA_PETSC_DIR}\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_PETSC_VERSION: ${UTOPIA_PETSC_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - slepc:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_SLEPC: ${UTOPIA_ENABLE_SLEPC}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_SLEPC_DIR: ${UTOPIA_SLEPC_DIR}\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_SLEPC_VERSION: ${UTOPIA_SLEPC_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - trilinos:\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_ENABLE_TRILINOS: ${UTOPIA_ENABLE_TRILINOS}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_TRILINOS_DIR: ${UTOPIA_TRILINOS_DIR}\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_TRILINOS_VERSION: ${UTOPIA_TRILINOS_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - blas:\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_ENABLE_BLAS: ${UTOPIA_ENABLE_BLAS}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_BLAS_DIR: ${UTOPIA_BLAS_DIR}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_BLAS_VERSION: ${UTOPIA_BLAS_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - lapack:\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_ENABLE_LAPACK: ${UTOPIA_ENABLE_LAPACK}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_LAPACK_DIR: ${UTOPIA_LAPACK_DIR}\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_LAPACK_VERSION: ${UTOPIA_LAPACK_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - umfpack:\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_ENABLE_UMFPACK: ${UTOPIA_ENABLE_UMFPACK}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_UMFPACK_DIR: ${UTOPIA_UMFPACK_DIR}\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_UMFPACK_VERSION: ${UTOPIA_UMFPACK_VERSION}\n")

  set(DEP_TABLE "${DEP_TABLE}  - yaml-cpp:\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_ENABLE_YAML_CPP: ${UTOPIA_ENABLE_YAML_CPP}\n")
  set(DEP_TABLE "${DEP_TABLE}    UTOPIA_YAML_CPP_DIR: ${UTOPIA_YAML_CPP_DIR}\n")
  set(DEP_TABLE
      "${DEP_TABLE}    UTOPIA_YAML_CPP_VERSION: ${UTOPIA_YAML_CPP_VERSION}\n")

  # OPTIONS

  set(DEP_TABLE "${DEP_TABLE}options:\n")
  getlistofvarsstartingwith("UTOPIA_" matchedVars)
  foreach(_var IN LISTS matchedVars)
    set(DEP_TABLE "${DEP_TABLE}  - ${_var}: ${${_var}}\n")
  endforeach()

  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfig.yaml" ${DEP_TABLE})
endmacro()
