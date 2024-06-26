file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/TestInstall)

add_custom_target(
    test_install
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${CMAKE_CURRENT_BINARY_DIR}/TestInstall/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove_directory
            ${CMAKE_CURRENT_BINARY_DIR}/TestInstall/CMakeFiles
    COMMAND
        ${CMAKE_COMMAND} -DUtopia_DIR=${CMAKE_INSTALL_PREFIX}/lib/cmake/
        ${CMAKE_SOURCE_DIR}/examples/usage_from_external_cmake_project
    COMMAND ${CMAKE_COMMAND} --build . --config $<IF:$<CONFIG:Debug>,Debug,Release>
    COMMAND ${CMAKE_CTEST_COMMAND} -V -C $<IF:$<CONFIG:Debug>,Debug,Release>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/TestInstall
    COMMENT "Testing the installation"
    VERBATIM)