# InstallYAMLCPP.cmake https://github.com/jbeder/yaml-cpp.git

include(ExternalProject)

if(UTOPIA_DEPENDENCIES_DIR)
  set(YAMLCPP_INSTALL_DIR ${UTOPIA_DEPENDENCIES_DIR}/yaml-cpp)
else()
  set(YAMLCPP_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../external/yaml-cpp)
endif()

set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
set(YAMLCPP_URL https://github.com/jbeder/yaml-cpp.git)
set(YAMLCPP_SOURCE_DIR "${STAGE_DIR}/yaml-cpp")

# ##############################################################################

ExternalProject_Add(
  yaml-cpp
  PREFIX "${STAGE_DIR}"
  GIT_REPOSITORY "${YAMLCPP_URL}"
  DOWNLOAD_DIR "${STAGE_DIR}"
  INSTALL_DIR "${YAMLCPP_INSTALL_DIR}"
  CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${YAMLCPP_INSTALL_DIR}"
  LOG_CONFIGURE 1
  LOG_BUILD 1
  BUILD_COMMAND ${CMAKE_COMMAND} -E echo "Starting $<CONFIG> build"
  COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config $<CONFIG>
  COMMAND ${CMAKE_COMMAND} -E echo "$<CONFIG> build complete")

set_target_properties(yaml-cpp PROPERTIES EXCLUDE_FROM_ALL TRUE)
