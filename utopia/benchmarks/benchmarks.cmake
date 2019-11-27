
set(UTOPIA_BENCH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)

list(APPEND BENCH_MODULES
    .
)

find_project_files(UTOPIA_BENCH_DIR BENCH_MODULES LOCAL_HEADERS LOCAL_SOURCES)
target_sources(utopia_bench PRIVATE ${LOCAL_SOURCES})
target_link_libraries(utopia_bench utopia)
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

# if(TARGET utopia_opencl)
#     target_link_libraries(utopia_bench utopia_opencl)
# endif()

# if(TARGET utopia_petsc)
#     target_link_libraries(utopia_bench utopia_petsc)
# endif()

# if(TARGET utopia_blas)
#     target_link_libraries(utopia_bench utopia_blas)
# endif()

# if(TARGET utopia_cuda_cxx)
#     target_link_libraries(utopia_bench utopia_cuda_cxx)
# endif()

# if(TARGET utopia_trilinos)
#     target_link_libraries(utopia_bench utopia_trilinos)
# endif()

# if(TARGET utopia_m3elinsol)
#     target_link_libraries(utopia_bench utopia_m3elinsol)
# endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")
