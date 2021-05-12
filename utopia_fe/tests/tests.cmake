# tests.cmake

if(UTOPIA_ENABLE_STK)
    list(APPEND TEST_MODULES stk)
endif()

if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND TEST_MODULES libmesh libmesh/old)
endif()

if(UTOPIA_ENABLE_MOONOLITH)
    list(APPEND TEST_MODULES moonolith)
endif()

if(UTOPIA_ENABLE_INTREPID2)
    list(APPEND TEST_MODULES intrepid2)
endif()

if(UTOPIA_ENABLE_MOONOLITH AND UTOPIA_ENABLE_STK)
    list(APPEND TEST_MODULES interop/moonolith_stk)
endif()

if(UTOPIA_ENABLE_INTREPID2 AND UTOPIA_ENABLE_STK)
    list(APPEND TEST_MODULES interop/stk_intrepid2)
endif()

if(UTOPIA_ENABLE_INTREPID2 AND UTOPIA_ENABLE_LIBMESH)
    list(APPEND TEST_MODULES interop/libmesh_intrepid2)
endif()

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")

find_project_files(${CMAKE_CURRENT_SOURCE_DIR}/tests "${TEST_MODULES}"
                   LOCAL_HEADERS LOCAL_SOURCES)

target_sources(utopia_fe_test PRIVATE ${LOCAL_SOURCES})
target_include_directories(utopia_fe_test
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests/generic)

foreach(MODULE ${TEST_MODULES})
    target_include_directories(
        utopia_fe_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests/${MODULE})
endforeach(MODULE)
