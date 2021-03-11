# tests.cmake

list(APPEND UTOPIA_FE_TEST_SOURCES
     ${CMAKE_SOURCE_DIR}/tests/utopia_NewLibMeshInterfaceTest.cpp)

target_sources(utopia_fe_test PRIVATE ${UTOPIA_FE_TEST_SOURCES})
