set(UTOPIA_BENCH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
set(UTOPIA_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

list(APPEND BENCH_MODULES .)


list(
    APPEND
    TEST_MODULES
    .
    test_problems
    test_problems/unconstrained_benchmark
    test_problems/constrained_benchmark
    test_problems/large_scale_benchmark
    test_problems/PTC_benchmark)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(${UTOPIA_BENCH_DIR} ${BENCH_MODULES} LOCAL_HEADERS
                   LOCAL_SOURCES)

target_sources(
  utopia_bench
  PRIVATE ${LOCAL_SOURCES}
  PRIVATE ${LOCAL_HEADERS})

foreach(MODULE ${TEST_MODULES})
  # deep dependency into the test module...
  target_include_directories(utopia_bench PRIVATE ${UTOPIA_BENCH_DIR}/${MODULE})
  target_include_directories(utopia_bench PRIVATE ${UTOPIA_TEST_DIR}/${MODULE})

endforeach()

if(UTOPIA_ENABLE_EIGEN_3)
  find_package(Eigen3)
  if(EIGEN3_FOUND)
    set(UTOPIA_ENABLE_EIGEN_3 ON)
    target_include_directories(utopia_bench PRIVATE ${EIGEN3_INCLUDE_DIR})
  endif()
endif()