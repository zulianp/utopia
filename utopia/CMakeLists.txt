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
include(cmake/UtopiaCompilerFeatures.cmake)
include(cmake/UtopiaGitVersion.cmake)
include(cmake/UtopiaClangFormat.cmake)
include(cmake/UtopiaFunctions.cmake)
include(cmake/CodeCoverage.cmake)
include(cmake/UtopiaDefaults.cmake)
include(cmake/UtopiaCompilerFeatures.cmake)
include(cmake/UtopiaDependencies.cmake)
include(cmake/UtopiaSubmodules.cmake)


if(UTOPIA_ENABLE_FLUYA_MODE)
  include(EnableFluyaMode)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

set(UTOPIA_ENABLE_CPP11 ON)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# LIBRARY

# Could take out backend. core modules
set(UTOPIA_MODULES interoperability solvers differentiable core ui utils)

foreach(MODULE ${UTOPIA_MODULES})
  add_subdirectory("${MODULE}")
endforeach(MODULE)

# Utopia Headers and Sources.
list(APPEND UTOPIA_HEADERS utopia.hpp utopia_Base.hpp utopia_Instance.hpp
     utopia_DeprecatedHeaders.hpp)

list(APPEND UTOPIA_SOURCES utopia.cpp utopia_Instance.cpp utopia_Base.cpp)

add_library(utopia ${UTOPIA_SOURCES} ${UTOPIA_HEADERS})

# We use this to force c++14
set_utopia_compiler_features(utopia)

utopia_link_default_targets(utopia)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# INCLUDE
foreach(include ${UTOPIA_BUILD_INCLUDES})
  target_include_directories(utopia BEFORE PUBLIC $<BUILD_INTERFACE:${include}>)
endforeach()

target_include_directories(
  utopia BEFORE
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:include>)

target_include_directories(utopia PUBLIC ${UTOPIA_DEP_INCLUDES})
message(STATUS "UTOPIA_DEP_INCLUDES: ${UTOPIA_DEP_INCLUDES}")

# ##############################################################################
# ##############################################################################
# ##############################################################################
# LINK

target_link_libraries(utopia PUBLIC ${UTOPIA_DEP_LIBRARIES})

# ##############################################################################
# ##############################################################################
# ##############################################################################
# INSTALL

# apps suite
add_executable(utopia_exec main.cpp)
target_link_libraries(utopia_exec PRIVATE utopia)
include(apps/apps.cmake)
utopia_link_default_targets(utopia_exec)

# message(STATUS "${UTOPIA_HEADERS}")

install(
  TARGETS utopia utopia_exec
  EXPORT UtopiaTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

if(UTOPIA_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

add_definitions(${UTOPIA_DEFS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")

# config is created based on the successful backend compilations
configure_file(${UTOPIA_ROOT_PATH}/utopia_Config.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_Config.hpp)
configure_file(${UTOPIA_ROOT_PATH}/utopia_Version.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_Version.hpp)

# write_compiler_detection_header( FILE utopia_compiler_features.hpp PREFIX
# UTOPIA COMPILERS GNU Clang MSVC Intel AppleClang FEATURES # cxx_final #already
# in c++11 cxx_noexcept #already in c++11 cxx_attribute_deprecated # in c++14
# ALLOW_UNKNOWN_COMPILERS)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia_Version.hpp
        DESTINATION include)

install(FILES ${UTOPIA_ROOT_PATH}/utopia.hpp DESTINATION include)
install(FILES ${UTOPIA_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/utopia_Config.hpp DESTINATION include)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# SUITES

# test suite
if(UTOPIA_ENABLE_TESTS)
  add_executable(utopia_test tests/test.cpp)
  utopia_link_default_targets(utopia_test)
  target_link_libraries(utopia_test PUBLIC utopia)
  include(tests/tests.cmake)

  enable_testing()
  add_custom_target(utest COMMAND time ${CMAKE_BINARY_DIR}/utopia_test -verbose)
  add_test(all_tests ${CMAKE_BINARY_DIR}/utopia_test -verbose)
  add_dependencies(utest utopia_test)
endif()

# bench suite Certain bench need headers from test
if(UTOPIA_ENABLE_BENCHMARK)
  add_executable(utopia_bench benchmarks/bench.cpp)
  utopia_link_default_targets(utopia_bench)
  target_link_libraries(utopia_bench PUBLIC utopia)
  include(benchmarks/benchmarks.cmake)

  add_custom_target(bench COMMAND time ${CMAKE_BINARY_DIR}/utopia_bench
                                  -verbose)
  add_dependencies(bench utopia_bench)
endif()

# # WHY? # Generate files useful for dependent build systems list(REVERSE
# UTOPIA_DEP_LIBRARIES) list(REMOVE_DUPLICATES UTOPIA_DEP_LIBRARIES)
# list(REVERSE UTOPIA_DEP_LIBRARIES)
set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")

# ##############################################################################
# ##############################################################################
# ##############################################################################
# MAKEFILE CONFIG

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
# ##############################################################################
# ##############################################################################
# ##############################################################################
# # MORE INSTALL and POST-INSTALL

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
              ${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfigVersion.cmake
        DESTINATION lib/cmake/)

include(cmake/PostInstall.cmake)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# EXTRA TARGETS

# ##############################################################################

if(UTOPIA_ENABLE_BENCHMARK AND UTOPIA_ENABLE_TESTS)
  add_custom_target(bench_and_test)
  add_dependencies(bench_and_test bench utest)
endif()

# ##############################################################################

if(UTOPIA_ENABLE_EXAMPLES
   AND UTOPIA_ENABLE_BENCHMARK
   AND UTOPIA_ENABLE_TESTS)
  add_custom_target(complete)
  add_dependencies(complete utopia utopia_test utopia_bench utopia_exec
                   utopia_examples)
endif()

print_dependency_table()
