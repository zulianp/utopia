utopia_add_library(utopia_polymorphic ${CMAKE_CURRENT_SOURCE_DIR} "." "utopia_core;utopia_utils;utopia_solvers;utopia_differentiable;utopia_ui;utopia_interoperability")

if(TRY_WITH_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(WITH_EIGEN_3 ON PARENT_SCOPE)
        target_include_directories(utopia_polymorphic PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

if(TARGET utopia_opencl)
    target_link_libraries(utopia_polymorphic PUBLIC utopia_opencl)
endif()

if(TARGET utopia_petsc)
    target_link_libraries(utopia_polymorphic PUBLIC utopia_petsc)
endif()

if(TARGET utopia_blas)
    target_link_libraries(utopia_polymorphic PUBLIC utopia_blas)
endif()

if(TARGET utopia_cuda_cxx)
    target_link_libraries(utopia_polymorphic PUBLIC utopia_cuda_cxx)
endif()

if(TARGET utopia_trilinos)
    target_link_libraries(utopia_polymorphic PUBLIC utopia_trilinos)
endif()

if(TARGET utopia_m3elinsol)
    target_link_libraries(utopia_polymorphic PUBLIC utopia_m3elinsol)
endif()
