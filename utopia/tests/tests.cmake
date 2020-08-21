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

# if(TARGET utopia_petsc) list(APPEND TEST_MODULES petsc) endif()

# if(TARGET utopia_trilinos) list(APPEND TEST_MODULES trilinos) endif()

# if(TARGET utopia_blas) list(APPEND TEST_MODULES blas) endif()

if(UTOPIA_PETSC)
    list(APPEND TEST_MODULES petsc)
endif()

if(UTOPIA_TRILINOS)
    list(APPEND TEST_MODULES trilinos)
endif()

if(UTOPIA_BLAS)
    list(APPEND TEST_MODULES blas)
endif()

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(${UTOPIA_TEST_DIR} "${TEST_MODULES}" LOCAL_HEADERS
                   LOCAL_SOURCES)

# message(STATUS "LOCAL_SOURCES=${LOCAL_SOURCES}
# LOCAL_HEADERS=${LOCAL_HEADERS}")

target_sources(utopia_test PRIVATE ${LOCAL_SOURCES})
target_link_libraries(utopia_test PRIVATE utopia)
utopia_link_default_targets(utopia_test)

target_include_directories(utopia_test PRIVATE ${UTOPIA_TEST_DIR})
target_include_directories(utopia_test PRIVATE .)
target_include_directories(utopia_test
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)

foreach(MODULE ${TEST_MODULES})
    target_include_directories(utopia_test PRIVATE ${UTOPIA_TEST_DIR}/${MODULE})
endforeach()

if(TRY_UTOPIA_WITH_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(UTOPIA_WITH_EIGEN_3 ON)
        # set(UTOPIA_WITH_EIGEN_3 ON PARENT_SCOPE)
        target_include_directories(utopia_test PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")
