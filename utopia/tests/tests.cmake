set(UTOPIA_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

list(
    APPEND
    TEST_MODULES
    .
    test_problems
    test_problems/unconstrained_benchmark
    test_problems/constrained_benchmark
    test_problems/large_scale_benchmark
    test_problems/PTC_benchmark)

list(APPEND TEST_MODULES deprecated)

if(UTOPIA_ENABLE_BENCHMARK)
    list(APPEND TEST_MODULES utest_bench)
endif()

#if(UTOPIA_ENABLE_PETSC)
#    list(APPEND TEST_MODULES petsc)
#endif()

if(UTOPIA_ENABLE_TRILINOS)
    list(APPEND TEST_MODULES trilinos)
endif()

if(UTOPIA_ENABLE_BLAS)
    list(APPEND TEST_MODULES blas)
endif()

if(UTOPIA_ENABLE_VC)
    list(APPEND TEST_MODULES vc)
endif()

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(${UTOPIA_TEST_DIR} "${TEST_MODULES}" LOCAL_HEADERS
                   LOCAL_SOURCES)

target_sources(utopia_test PRIVATE ${LOCAL_SOURCES} PRIVATE ${LOCAL_HEADERS})

target_include_directories(utopia_test
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)

foreach(MODULE ${TEST_MODULES})
    target_include_directories(utopia_test PRIVATE ${UTOPIA_TEST_DIR}/${MODULE})
endforeach()

if(UTOPIA_ENABLE_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(UTOPIA_ENABLE_EIGEN_3 ON)
        target_include_directories(utopia_test PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()
