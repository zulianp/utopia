# OLD MAIN CMAKE FILE FOR REFERENCE


cmake_minimum_required(VERSION 3.13...3.25)
cmake_policy(SET CMP0079 NEW)
include(cmake/UtopiaOptions.cmake)
project(
  utopia
  VERSION ${UTOPIA_VERSION}
  LANGUAGES "CXX" "C"
  HOMEPAGE_URL "https://bitbucket.org/zulianp/utopia/"
  DESCRIPTION
    "Utopia is an open-source C++ library for parallel non-linear multilevel solution strategies"
)

# ##############################################################################
# ######################## UTOPIA DEFAULT INCLUDES
# ##############################################################################
# ##############################################################################

include(CMakePackageConfigHelpers)
include(cmake/UtopiaGitVersion.cmake)
include(cmake/UtopiaClangFormat.cmake)
include(cmake/UtopiaFunctions.cmake)
include(cmake/CodeCoverage.cmake)
include(cmake/UtopiaDefaults.cmake)

if(UTOPIA_ENABLE_FLUYA_MODE)
  include(EnableFluyaMode)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

set(UTOPIA_WITH_CPP11 ON)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# LIBRARY



# Utopia Headers and Sources.
list(APPEND UTOPIA_HEADERS utopia.hpp utopia_Base.hpp utopia_Instance.hpp
     utopia_DeprecatedHeaders.hpp)

list(APPEND UTOPIA_SOURCES utopia.cpp utopia_Instance.cpp utopia_Base.cpp)

add_library(utopia ${UTOPIA_SOURCES})


# TEMPORAY PLACEMENT OF FUNCTION.
function(utopia_add_library libraryRootDir subDirs)
  set(THIS_MODULE_SOURCES "")
  set(THIS_MODULE_HEADERS "")
  find_project_files(${libraryRootDir} "${subDirs}" THIS_MODULE_HEADERS
                     THIS_MODULE_SOURCES)
  target_sources(
    utopia
    PRIVATE ${THIS_MODULE_SOURCES}
    PRIVATE ${THIS_MODULE_HEADERS})

  install(FILES ${THIS_MODULE_HEADERS} DESTINATION include)
  foreach(MODULE ${subDirs})
    target_include_directories(
      utopia BEFORE PUBLIC $<BUILD_INTERFACE:${libraryRootDir}/${MODULE}>)
  endforeach(MODULE)
endfunction()

# ALSO TEMPORARY
include(cmake/UtopiaDependencies.cmake)

# core modules
set(UTOPIA_MODULES
    interoperability
    backend
    solvers
    differentiable
    core
    ui
    utils)

foreach(MODULE ${UTOPIA_MODULES})
  add_subdirectory(${MODULE})
endforeach(MODULE)

# if(UTOPIA_ENABLE_KOKKOS_SIMD) set(UTOPIA_WITH_KOKKOS_SIMD TRUE)

# target_include_directories(utopia PUBLIC
# $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/simd-math>) endif()

# We use this to force c++14
set_utopia_compiler_features(utopia)

utopia_link_default_targets(utopia)
target_include_directories(
  utopia BEFORE PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
                       $<INSTALL_INTERFACE:include>)
target_include_directories(
  utopia BEFORE PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(
  TARGETS utopia
  EXPORT UtopiaTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

add_subdirectory(examples)
add_subdirectory(backend/polymorphic)


add_definitions(${UTOPIA_DEFS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")

# config is created based on the successful backend compilations
configure_file(${UTOPIA_ROOT_PATH}/utopia_Config.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_Config.hpp)
configure_file(${UTOPIA_ROOT_PATH}/utopia_Version.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_Version.hpp)

# ##############################################################################

# Deprecated since version 3.20 include(WriteCompilerDetectionHeader)

write_compiler_detection_header( FILE utopia_compiler_features.hpp PREFIX
UTOPIA COMPILERS GNU Clang MSVC Intel AppleClang FEATURES # cxx_final #already
in c++11 cxx_noexcept #already in c++11 cxx_attribute_deprecated # in c++14
ALLOW_UNKNOWN_COMPILERS)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia_compiler_features.hpp
              ${CMAKE_CURRENT_BINARY_DIR}/utopia_Version.hpp
        DESTINATION include)

install(FILES ${UTOPIA_ROOT_PATH}/utopia.hpp DESTINATION include)
install(FILES ${UTOPIA_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/utopia_Config.hpp DESTINATION include)

# ##############################################################################

# test suite
add_executable(utopia_test EXCLUDE_FROM_ALL tests/test.cpp)
utopia_link_default_targets(utopia_test)
target_link_libraries(utopia_test PRIVATE utopia)
include(tests/tests.cmake)

# bench suite
add_executable(utopia_bench EXCLUDE_FROM_ALL benchmarks/bench.cpp)
utopia_link_default_targets(utopia_bench)
target_link_libraries(utopia_bench PRIVATE utopia)
include(benchmarks/benchmarks.cmake)

# apps suite
add_executable(utopia_exec EXCLUDE_FROM_ALL main.cpp)
utopia_link_default_targets(utopia_exec)
target_link_libraries(utopia_exec PRIVATE utopia)
include(apps/apps.cmake)


# Generate files useful for dependent build systems
list(REVERSE UTOPIA_DEP_LIBRARIES)
list(REMOVE_DUPLICATES UTOPIA_DEP_LIBRARIES)
list(REVERSE UTOPIA_DEP_LIBRARIES)
set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")


 # Make file stuff 
append_list_to_string_for_makefile(
  CONFIG_CONTENT "UTOPIA_INCLUDES=-I${CMAKE_INSTALL_PREFIX}/include " "\n" "-I"
  "${UTOPIA_BUILD_INCLUDES}")
append_list_to_string_for_makefile(
  CONFIG_CONTENT "UTOPIA_LIBRARIES=-L${CMAKE_INSTALL_PREFIX}/lib " " -lutopia"
  "-l" "${UTOPIA_COMPONENTS}")
append_list_to_string_for_makefile(CONFIG_CONTENT " " "\n" ""
                                   "${UTOPIA_DEP_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_DEFS=" "\n" ""
                                   "${UTOPIA_DEFS}")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_CXX_FLAGS=${CMAKE_CXX_FLAGS}\n")
set(CONFIG_CONTENT
    "${CONFIG_CONTENT}UTOPIA_CXX_COMPILER=${CMAKE_CXX_COMPILER}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_C_COMPILER=${CMAKE_C_COMPILER}\n")

file(WRITE ${CMAKE_BINARY_DIR}/utopia-config.makefile ${CONFIG_CONTENT})

set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(
  CONFIG_CONTENT
  "LIST(APPEND UTOPIA_INCLUDES\n\t${CMAKE_INSTALL_PREFIX}/include\n" "\n)\n\n"
  "${UTOPIA_BUILD_INCLUDES}")
append_list_to_string_for_makefile(
  CONFIG_CONTENT
  "LIST(APPEND UTOPIA_LIBRARIES\n\t-L${CMAKE_INSTALL_PREFIX}/lib\n\t-lutopia"
  "" "\n\t-l" "${UTOPIA_COMPONENTS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "\n" ")\n\n"
                                "${UTOPIA_DEP_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_DEFS\n"
                                "\n)\n\n" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_FLAGS\n"
                                "\n)\n\n" "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(
  CONFIG_CONTENT "LIST(APPEND UTOPIA_C_COMPILER\n" "\n)\n\n"
  "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(
  CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_COMPILER\n" "\n)\n\n"
  "${CMAKE_CXX_COMPILER}")
file(WRITE ${CMAKE_BINARY_DIR}/utopia-config.cmake ${CONFIG_CONTENT})
list(JOIN CONFIG_CONTENT "\n" OLD_CMAKE_CONFIG_CONTENT)

install(FILES ${CMAKE_BINARY_DIR}/utopia-config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia-config.cmake DESTINATION config)

# new style exported targets, this allows to use target_link_libraries(my_target
# Utopia::utopia)
install(
  EXPORT UtopiaTargets
  FILE UtopiaTargets.cmake
  NAMESPACE Utopia::
  DESTINATION lib/cmake)

# generate the config file that is includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/UtopiaConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfig.cmake"
  INSTALL_DESTINATION "lib/cmake")
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfigVersion.cmake"
  VERSION "${Utopia_VERSION_MAJOR}.${Utopia_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)
# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfig.cmake
        DESTINATION lib/cmake/)

get_target_property(XX_LINK_DIRECTORIES utopia_exec LINK_DIRECTORIES)
get_target_property(XX_LINK_LIBRARIES utopia_exec LINK_INTERFACE_LIBRARIES)
get_target_property(XX_INCLUDE_DIRECTORIES utopia_exec
                    INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "LINK_DIRECTORIES=${XX_LINK_DIRECTORIES}")
message(STATUS "LINK_LIBRARIES=${XX_LINK_LIBRARIES}")
message(STATUS "INCLUDE_DIRECTORIES=${XX_INCLUDE_DIRECTORIES}")

print_dependency_table()

# find_package(OpenMP) if(OPENMP_FOUND) SET (CMAKE_Fortran_FLAGS
# "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}") SET (CMAKE_CXX_FLAGS
# "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}") SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}
# ${OpenMP_C_FLAGS}") set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS}
# ${OpenMP_CXX_FLAGS}") endif()

# ##############################################################################

# enable_testing()
add_custom_target(utest COMMAND time ${CMAKE_BINARY_DIR}/utopia_test -verbose)
# add_test(all_tests ${CMAKE_BINARY_DIR}/utopia_test -verbose)
add_dependencies(utest utopia_test)

# ##############################################################################

add_custom_target(bench COMMAND time ${CMAKE_BINARY_DIR}/utopia_bench -verbose)
add_dependencies(bench utopia_bench)

# ##############################################################################

add_custom_target(bench_and_test)
add_dependencies(bench_and_test bench utest)

# ##############################################################################

add_custom_target(complete)
add_dependencies(complete utopia utopia_test utopia_bench utopia_exec
                 utopia_examples)

# ##############################################################################
target_format_setup("utopia;utopia_test;utopia_bench;utopia_exec")
# ##############################################################################

# During the cmake step some test programs can crash and leave a core-dump file
# behind This file has unfortunately the same name as the module `core`.
# Therefore at the end when CMake tries to generate the output directories
# (which happens after the whole CMakeLists.txt has been parsed), we must make
# sure that there is no such file with the name `core`. One program that is
# known to crash on Piz Daint is the PETSc test run, because MPI cannot be
# initialised when it is not running under an srun environment (i.e. you are
# compiling on a login node).
if(EXISTS ${CMAKE_BINARY_DIR}/core AND NOT IS_DIRECTORY
                                       ${CMAKE_BINARY_DIR}/core)
  file(REMOVE ${CMAKE_BINARY_DIR}/core)
endif()
