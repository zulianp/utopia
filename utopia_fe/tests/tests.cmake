# tests.cmake

# Added intrepid2 since utopia_stk_restartTest need intrepid.
if(UTOPIA_ENABLE_STK AND UTOPIA_ENABLE_INTREPID2)
    list(APPEND TEST_MODULES stk)
endif()

if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND TEST_MODULES libmesh)

    # if(UTOPIA_ENABLE_LIBMESH_DEPRECATED)
    #     list(APPEND TEST_MODULES libmesh/deprecated)
    # endif()
endif()

if(UTOPIA_ENABLE_MOONOLITH)
    list(APPEND TEST_MODULES moonolith)
endif()

if(UTOPIA_ENABLE_INTREPID2)
    list(APPEND TEST_MODULES intrepid2)
endif()

if(UTOPIA_ENABLE_MARS)
    list(APPEND TEST_MODULES mars)
endif()

if(UTOPIA_ENABLE_PETSC)
    list(APPEND TEST_MODULES petsc)
endif()

if(UTOPIA_ENABLE_MOONOLITH AND UTOPIA_ENABLE_STK)
    list(APPEND TEST_MODULES interop/moonolith_stk)
endif()

if(UTOPIA_ENABLE_INTREPID2 AND UTOPIA_ENABLE_STK)
    list(APPEND TEST_MODULES interop/stk_intrepid2)
endif()

if(UTOPIA_ENABLE_INTREPID2 AND UTOPIA_ENABLE_PETSC)
    list(APPEND TEST_MODULES interop/petsc_intrepid2)
endif()

if(UTOPIA_ENABLE_INTREPID2 AND UTOPIA_ENABLE_ARBORX)
    list(APPEND TEST_MODULES interop/intrepid2_arborx)
endif()

if(UTOPIA_ENABLE_INTREPID2 AND UTOPIA_ENABLE_LIBMESH AND UTOPIA_ENABLE_LIBMESH_KOKKOS)
    list(APPEND TEST_MODULES interop/libmesh_kokkos)
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
