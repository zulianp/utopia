utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")

if(UTOPIA_ENABLE_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(UTOPIA_WITH_EIGEN_3 ON PARENT_SCOPE)
        target_include_directories(utopia PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

if(TARGET utopia_opencl)
    target_link_libraries(utopia PUBLIC utopia_opencl)
endif()

if(TARGET utopia_petsc)
    target_link_libraries(utopia PUBLIC utopia_petsc)
endif()

if(TARGET utopia_blas)
    target_link_libraries(utopia PUBLIC utopia_blas)
endif()

if(TARGET utopia_cuda_cxx)
    target_link_libraries(utopia PUBLIC utopia_cuda_cxx)
endif()

if(TARGET utopia_trilinos)
    target_link_libraries(utopia PUBLIC utopia_trilinos)
endif()

if(TARGET utopia_m3elinsol)
    target_link_libraries(utopia PUBLIC utopia_m3elinsol)
endif()
