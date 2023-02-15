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
include(cmake/UtopiaDependencies.cmake)

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
  add_subdirectory("${MODULE}")
endforeach(MODULE)

# Utopia Headers and Sources.
list(APPEND UTOPIA_HEADERS utopia.hpp utopia_Base.hpp utopia_Instance.hpp
     utopia_DeprecatedHeaders.hpp)

list(APPEND UTOPIA_SOURCES utopia.cpp utopia_Instance.cpp utopia_Base.cpp)

add_library(utopia ${UTOPIA_SOURCES} ${UTOPIA_HEADERS})

# PRINT OUT HEADERS message(STATUS "UTOPIA_HEADERS") foreach(header
# ${UTOPIA_HEADERS}) message(STATUS "${header}") endforeach(header)

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
message("UTOPIA_DEP_INCLUDES: ${UTOPIA_DEP_INCLUDES}")
message("UTOPIA_DEP_LIBRARIES: ${UTOPIA_DEP_LIBRARIES}")

target_link_libraries(utopia PUBLIC ${UTOPIA_DEP_LIBRARIES})

# ##############################################################################
# ##############################################################################
# ##############################################################################
# INSTALL

install(
  TARGETS utopia
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



# ##############################################################################
# ##############################################################################
# ##############################################################################
# EXTRA TARGETS


