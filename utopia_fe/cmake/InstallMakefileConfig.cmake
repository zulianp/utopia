# InstallMakefileConfig.cmake

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/InstallMakefileConfig)

add_custom_target(
  install_makefile_config
  COMMAND ${CMAKE_COMMAND} -E remove -f
          ${CMAKE_CURRENT_BINARY_DIR}/InstallMakefileConfig/CMakeCache.txt
  COMMAND ${CMAKE_COMMAND} -E remove_directory
          ${CMAKE_CURRENT_BINARY_DIR}/InstallMakefileConfig/CMakeFiles
  COMMAND
    ${CMAKE_COMMAND} -DUtopiaFE_DIR=${CMAKE_INSTALL_PREFIX}/lib/cmake/
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    ${CMAKE_SOURCE_DIR}/cmake/utils
  COMMAND ${CMAKE_COMMAND} --build . --target all
  COMMAND ${CMAKE_COMMAND} -P cmake_install.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/InstallMakefileConfig
  COMMENT "Installing configuration for makefile users."
  VERBATIM)
