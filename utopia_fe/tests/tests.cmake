# tests.cmake

if(UTOPIA_ENABLE_STK)
    # if(TARGET utopia_stk)
    list(APPEND TEST_MODULES stk)
    target_link_libraries(utopia_fe_test utopia_stk)
    # endif()
endif()

if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND TEST_MODULES libmesh)
endif()

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")

find_project_files(${CMAKE_CURRENT_SOURCE_DIR}/tests "${TEST_MODULES}"
                   LOCAL_HEADERS LOCAL_SOURCES)

target_sources(utopia_fe_test PRIVATE ${LOCAL_SOURCES})
