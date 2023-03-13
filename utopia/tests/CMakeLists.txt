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

if(UTOPIA_ENABLE_PETSC)
    list(APPEND TEST_MODULES petsc)
endif()

if(UTOPIA_ENABLE_TRILINOS)
    list(APPEND TEST_MODULES trilinos)
endif()

if(UTOPIA_ENABLE_BLAS)
    list(APPEND TEST_MODULES blas)
endif()

if(UTOPIA_ENABLE_VC)
    list(APPEND TEST_MODULES vc)
endif()

scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${TEST_MODULES}"
                 UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)