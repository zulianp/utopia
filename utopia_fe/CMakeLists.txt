cmake_minimum_required(VERSION 3.13...3.25)
cmake_policy(SET CMP0079 NEW)
project(
  utopia_fe
  VERSION 0.1.0
  LANGUAGES "CXX" "C"
  HOMEPAGE_URL "https://bitbucket.org/zulianp/utopia"
  DESCRIPTION "Library for portable finite element discretizations.")

# ##############################################################################
# ##############################################################################
# ##############################################################################
# MODES

if(UTOPIA_ENABLE_FLUYA_MODE)
  include(cmake/EnableFluyaMode.cmake)
endif()

if(UTOPIA_ENABLE_FRANETG_MODE)
  include(cmake/EnableFranetgMode.cmake)
endif()

if(UTOPIA_ENABLE_AVFLOW_MODE)
  include(cmake/EnableAvflowMode.cmake)
endif()
# ##############################################################################
# ######################## UTOPIA_FE DEFAULT INCLUDES
# ##############################################################################
# ##############################################################################

include(cmake/CodeCoverage.cmake)
include(cmake/UtopiaFEOptions.cmake)
include(cmake/UtopiaFEFunctions.cmake)
include(cmake/UtopiaFEDefaults.cmake)
include(cmake/UtopiaFEDependencies.cmake)
include(cmake/UtopiaCompilerFeatures.cmake)

# file(GLOB_RECURSE FUNC_LIB "{UTOPIA_FE_ROOT_PATH}/cmake/utils")
# foreach(funcLib ${FUNC_LIB}) include(${funcLib}) endforeach(funcLib)

# ##############################################################################

list(APPEND UTOPIA_FE_MODULES frontend interoperability base core)

if(UTOPIA_ENABLE_LIBMESH)
  list(APPEND UTOPIA_FE_MODULES backend/libmesh edsl)
endif()

if(UTOPIA_ENABLE_MOONOLITH)
  list(APPEND UTOPIA_FE_MODULES backend/moonolith)
endif()

if(UTOPIA_ENABLE_STK)
  list(APPEND UTOPIA_FE_MODULES backend/stk)
endif()

if(UTOPIA_ENABLE_PETSCDM)
  list(APPEND UTOPIA_FE_MODULES backend/petsc)
endif()

if(UTOPIA_ENABLE_SFEM)
  list(APPEND UTOPIA_FE_MODULES backend/sfem)
endif()

# #
# ##############################################################################

if(UTOPIA_ENABLE_INTREPID2)
  list(APPEND UTOPIA_FE_MODULES backend/intrepid2)
endif()

if(UTOPIA_ENABLE_INTREPID2
   OR UTOPIA_ENABLE_MARS
   OR UTOPIA_ENABLE_LIBMESH_KOKKOS)
  list(APPEND UTOPIA_FE_MODULES backend/kokkos)
endif()

if(UTOPIA_ENABLE_MARS)
  list(APPEND UTOPIA_FE_MODULES backend/mars)
endif()

# if(UTOPIA_ENABLE_ARBORX) list(APPEND UTOPIA_FE_MODULES backend/arborx) endif()

foreach(MODULE ${UTOPIA_FE_MODULES})
  add_subdirectory("${MODULE}")
endforeach(MODULE)

# ##############################################################################

list(APPEND UTOPIA_FE_SOURCES utopia_fe.cpp)
list(APPEND UTOPIA_FE_HEADERS utopia_fe.hpp)

add_library(utopia_fe ${UTOPIA_FE_SOURCES} ${UTOPIA_FE_HEADERS})

if(UTOPIA_FE_TARGET_DEPENDENCIES)
  add_dependencies(utopia_fe ${UTOPIA_FE_TARGET_DEPENDENCIES})
endif()

if(UTOPIA_FE_SUBMODULES)
  target_link_libraries(utopia_fe PUBLIC ${UTOPIA_FE_SUBMODULES})
endif()

set_utopia_compiler_features(utopia_fe)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# INCLUDE

foreach(include ${UTOPIA_FE_BUILD_INCLUDES})
  target_include_directories(utopia_fe BEFORE
                             PUBLIC $<BUILD_INTERFACE:${include}>)
endforeach()

target_include_directories(
  utopia_fe BEFORE
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:include>)

list(REMOVE_DUPLICATES UTOPIA_FE_DEP_INCLUDES)
list(REMOVE_DUPLICATES UTOPIA_FE_DEP_LIBRARIES)

target_include_directories(utopia_fe PUBLIC ${UTOPIA_FE_DEP_INCLUDES})

# ##############################################################################
# ##############################################################################
# ##############################################################################
# LINK

target_link_libraries(utopia_fe PUBLIC ${UTOPIA_FE_DEP_LIBRARIES})

# ##############################################################################
# INSTALL

add_executable(utopia_fe_exec main.cpp)
target_link_libraries(utopia_fe_exec PRIVATE utopia_fe)
include(apps/apps.cmake)
utopia_fe_link_default_targets(utopia_fe_exec)

install(
  TARGETS utopia_fe utopia_fe_exec
  EXPORT UtopiaFETargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(FILES utopia_fe.hpp DESTINATION include)
configure_file(${UTOPIA_FE_ROOT_PATH}/utopia_fe_config.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp)
install(FILES ${UTOPIA_FE_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp DESTINATION include)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# SUITES

# test suite
if(UTOPIA_ENABLE_TESTS)
  add_executable(utopia_fe_test test.cpp)
  utopia_fe_link_default_targets(utopia_fe_test)
  target_link_libraries(utopia_fe_test PRIVATE utopia_fe)
  include(tests/tests.cmake)
  # target_link_libraries(utopia_fe_exec PRIVATE utopia_fe ${TEST_MODULES})
endif()

# # set(TOP_LEVEL 1)

# #
# --------------------------------------------------------------------------------------------------
# #
# --------------------------------------------------------------------------------------------------

# if(UTOPIA_FE_LIBRARIES_TRAILING) list(APPEND UTOPIA_FE_LIBRARIES
# ${UTOPIA_FE_LIBRARIES_TRAILING}) endif()

# Generate files useful for dependent build systems set(CONFIG_CONTENT
# "#Automatically generated by utopia CMakeLists.txt\n#Includes the content of
# utopia-config.makefile generated by utopia\n\n" )

# MAKEFILE CONFIG

# list(APPEND UTOPIA_FE_INCLUDES ${UTOPIA_INCLUDES})

append_list_to_string_for_makefile(
  CONFIG_CONTENT
  "UTOPIA_FE_INCLUDES=-I${CMAKE_INSTALL_PREFIX}/include -I${LIBMESH_INCLUDE_DIR}
" " " "-I" "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" "\n" "-I"
                                   "${MOONOLITH_INCLUDES}")

# message("MOONOLITH_INCLUDES=${MOONOLITH_INCLUDES}")

append_list_to_string_for_makefile(
  CONFIG_CONTENT
  "UTOPIA_FE_LIBRARIES=-L${CMAKE_INSTALL_PREFIX}/lib -lutopia_fe " " " "-l"
  "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " ""
                                   "${LIBMESH_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " ""
                                   "${UTOPIA_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " ""
                                   "${MOONOLITH_LIBRARIES}")

set(CONFIG_CONTENT "${CONFIG_CONTENT}\n")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_DEFS=" "\n" ""
                                   "${UTOPIA_DEFS}")

set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_CXX_FLAGS=${UTOPIA_CXX_FLAGS}\n")
set(CONFIG_CONTENT
    "${CONFIG_CONTENT}UTOPIA_FE_CXX_COMPILER=${UTOPIA_CXX_COMPILER}\n")
set(CONFIG_CONTENT
    "${CONFIG_CONTENT}UTOPIA_FE_C_COMPILER=${UTOPIA_C_COMPILER}\n")

file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile ${CONFIG_CONTENT})

set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(
  CONFIG_CONTENT "LIST(APPEND
UTOPIA_FE_INCLUDES\n\t${CMAKE_INSTALL_PREFIX}/include\n" ""
  "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" "" "${LIBMESH_INCLUDE_DIR}")

if(MOONOLITH_INCLUDES)
  append_list_to_string_for_cmake(CONFIG_CONTENT "" ")\n\n"
                                  ${MOONOLITH_INCLUDES})
else()
  set(CONFIG_CONTENT "${CONFIG_CONTENT})")
endif()

append_list_to_string_for_makefile_with_postfix(
  CONFIG_CONTENT "LIST(APPEND
UTOPIA_FE_LIBRARIES\n\t-L${CMAKE_INSTALL_PREFIX}/lib\n\t-lutopia_fe" ""
  "\n\t\"-l" "\"" "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "" "" "\n\t\""
                                                "\"" "${MOONOLITH_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "" "" "\n\t\""
                                                "\"" "${LIBMESH_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "\n" "\n)\n\n"
                                "${UTOPIA_LIBRARIES}")

append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_DEFS\n"
                                "" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" "\n)\n\n" "${UTOPIA_FE_DEFS}")
append_list_to_string_for_cmake(
  CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_CXX_FLAGS\n" "\n)\n\n"
  "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND
UTOPIA_FE_C_COMPILER\n" "\n)\n\n" "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND
UTOPIA_FE_CXX_COMPILER\n" "\n)\n\n" "${CMAKE_CXX_COMPILER}")

file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake ${CONFIG_CONTENT})
list(JOIN CONFIG_CONTENT "\n" OLD_CMAKE_CONFIG_CONTENT)

# ##############################################################################
# ##############################################################################
# ##############################################################################
# # MORE INSTALL and POST-INSTALL

install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake DESTINATION config)

# new style exported targets, this allows to use target_link_libraries(my_target
# UtopiaFE::utopia_fe)

install(
  EXPORT UtopiaFETargets
  FILE UtopiaFETargets.cmake
  NAMESPACE UtopiaFE::
  DESTINATION lib/cmake)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/UtopiaFEConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/UtopiaFEConfig.cmake"
  INSTALL_DESTINATION "lib/cmake")
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/UtopiaFEConfigVersion.cmake"
  VERSION "${UtopiaFE_VERSION_MAJOR}.${UtopiaFE_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)

# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UtopiaFEConfig.cmake
        DESTINATION lib/cmake/)

if(UTOPIA_ENABLE_LIBMESH AND UTOPIA_ENABLE_LIBMESH_DEPRECATED)
  include(deprecated/apps/apps.cmake)
endif()

# POST-INSTALL
# ##############################################################################
include(cmake/PostInstall.cmake)

# EXTRA TARGETS
# ##############################################################################

print_dependency_table()
log_dependency_table()
