# option(UTOPIA_ENABLE_VC "SIMD wrapper" ON)
option(UTOPIA_ENABLE_VC "SIMD wrapper" OFF)

option(UTOPIA_STATIC_DEPENDENCIES_ONLY
       "Allows to restrict the linking to static libraries" OFF)

option(UTOPIA_ENABLE_CODE_COVERAGE "Enable coverage reporting" OFF)

option(UTOPIA_ENABLE_ENV_READ "Enable utopia to look at enviroment variables" ON)   

option(UTOPIA_ENABLE_EXAMPLES "Enable utopia examples." OFF)
option(UTOPIA_ENABLE_TESTS "Enable utopia tests." ON)
option(UTOPIA_ENABLE_BENCHMARK "Enable utopia benchmarks." ON)
option(UTOPIA_ENABLE_MPI "Enable the cuda backend" ON)
# option(UTOPIA_ENABLE_OPENCL "Enable the opencl backend" OFF)
option(UTOPIA_ENABLE_PETSC "Enable the petsc backend" ON)
option(UTOPIA_ENABLE_SLEPC "Enable the slepc backend" OFF)
option(UTOPIA_ENABLE_METIS "Enable the Metis backend" OFF)
option(UTOPIA_ENABLE_POLYMORPHIC "Enable polymorphic" ON)
option(UTOPIA_ENABLE_PARMETIS "Enable the ParMetis backend" OFF)
option(UTOPIA_ENABLE_TRILINOS "Enable the Trilinos backend" ON)
# option(UTOPIA_ENABLE_KOKKOS_SIMD "Enable kokkos intriniscs wrapper" OFF)
option(UTOPIA_ENABLE_BLAS "Enable the blas backend" ON)
option(UTOPIA_ENABLE_MARS "Enable the mars backend" OFF)

option(UTOPIA_ENABLE_TRACE "enables utopia tracing facilities for regions" OFF)
option(UTOPIA_ENABLE_TRACE_EXPR "enables utopia tracing facilities for every
expression" OFF)

option(UTOPIA_ENABLE_LOCAL_DEPENDENCIES_INSTALL "Sets utopia to use current project installs of petsc and trilinos" OFF)

option(UTOPIA_ENABLE_NO_ALLOC_REGIONS"enables utopia alloc-check facilities"
       OFF)
option(UTOPIA_ENABLE_SANITIZER "Enable clang -fsanitize=address flag" OFF)
# option(UTOPIA_ENABLE_PASSO_EXTENSIONS "Enable non-standard petsc solvers
# developed in the PASSO library" OFF)
option(UTOPIA_ENABLE_CXX14_FEATURES "Enable usage of cxx14 standard" ON)
option(UTOPIA_ENABLE_CXX17_FEATURES "Enable usage of cxx17 standard" ON)

option(UTOPIA_ENABLE_MOOSE_ENV_MODE
       "Allows to use the moose installation compilers" OFF)
option(UTOPIA_ENABLE_DEPRECATED_API
       "Decprecated functionality that will be removed in the future" ON)
option(UTOPIA_ENABLE_TINY_EXPR "String expressions support" ON)
option(UTOPIA_ENABLE_RAPIDJSON "Enable support for JSON input files" ON)
option(UTOPIA_ENABLE_YAML_CPP "Enable YAML support" ON)

option(UTOPIA_ENABLE_CLANG_TIDY "Use clang tidy static analyzer" OFF)
option(
  UTOPIA_ENABLE_PETSC_DM_PLEX
  "For the moment is to allow the \"make petsc\" command to install DMPlex dependencies"
  OFF)
option(UTOPIA_PULL_REQUEST_MODE
       "Mode for ensuring that pull request standards are respected" OFF)

option(UTOPIA_ENABLE_EIGEN_3 "Look for eigen for comparing performance" OFF)
option(UTOPIA_BUILD_DOCUMENTATION
       "Use Doxygen to create the HTML based API documentation" ON)

option(UTOPIA_ENABLE_SCRIPTING
       "Enable exports for other languages (C, Python, etc...)" OFF)

option(UTOPIA_INSTALL_TRILINOS "Install trilinos directly" ON)
option(UTOPIA_INSTALL_PETSC "Install petsc directly" ON)
option(UTOPIA_INSTALL_PETSC_DEBUG "Install petsc directly" OFF)
option(
  UTOPIA_ENABLE_GPERFTOOLS
  "use Google perf tools (https://gperftools.github.io/gperftools/cpuprofile.html)"
  OFF)

option(UTOPIA_SPACK_MODE "Spack env hacks" OFF)

option(BUILD_SHARED_LIBS "Build shared libraries instead of static" OFF)

# Setting this to be on, otherwise will not compile since Kokkos_DefaultNode.hpp
# does not exist anymore.
option(UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE
       "Removed trilinos deprecated code" ON)

option(UTOPIA_ENABLE_FLUYA_MODE "Create utopia configuration required by Fluya"
       OFF)

option(UTOPIA_ENABLE_ISOLVER "Enable ISolver" OFF)

set(UTOPIA_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_ROOT_PATH}/cmake")

set(CMAKE_C_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
          FORCE)

set(CMAKE_CXX_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING
          "Flags used by the C++ compiler during AddressSanitizer builds."
          FORCE)

# ##############################################################################
# ######################## UTOPIA VERSION INFORMATION
# ##############################################################################
# ##############################################################################

set(UTOPIA_VERSION_MAJOR 0)
set(UTOPIA_VERSION_MINOR 1)
set(UTOPIA_VERSION_PATCH 0)
set(UTOPIA_VERSION
    "${UTOPIA_VERSION_MAJOR}.${UTOPIA_VERSION_MINOR}.${UTOPIA_VERSION_PATCH}")
# ##############################################################################

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

  message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

# To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)

set(CMAKE_CXX_EXTENSIONS OFF)
