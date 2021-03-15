# tests.cmake

if(UTOPIA_ENABLE_STK)
    list(APPEND TEST_MODULES stk)
endif()

if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND TEST_MODULES libmesh)
endif()

if(UTOPIA_ENABLE_MOONOLITH)
    list(APPEND TEST_MODULES moonolith)
endif()

if(UTOPIA_ENABLE_MOONOLITH AND UTOPIA_ENABLE_STK)
    list(APPEND TEST_MODULES interop)
endif()

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")

find_project_files(${CMAKE_CURRENT_SOURCE_DIR}/tests "${TEST_MODULES}"
                   LOCAL_HEADERS LOCAL_SOURCES)

target_sources(utopia_fe_test PRIVATE ${LOCAL_SOURCES})
