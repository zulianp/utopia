set(UTOPIA_BENCH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)

list(APPEND BENCH_MODULES
    .
)



set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(UTOPIA_BENCH_DIR BENCH_MODULES LOCAL_HEADERS LOCAL_SOURCES)
target_sources(utopia_bench PRIVATE ${LOCAL_SOURCES})
target_link_libraries(utopia_bench utopia)
utopia_link_default_targets(utopia_bench)

target_include_directories(utopia_bench PRIVATE ${UTOPIA_BENCH_DIR})
target_include_directories(utopia_bench PRIVATE .)
target_include_directories(utopia_bench PRIVATE ${BENCH_MODULES})

if(TRY_WITH_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(WITH_EIGEN_3 ON)
        set(WITH_EIGEN_3 ON PARENT_SCOPE)
        target_include_directories(utopia_bench PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")
