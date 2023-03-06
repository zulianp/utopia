# utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")
scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "." UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)

if(UTOPIA_ENABLE_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(UTOPIA_ENABLE_EIGEN_3 ON PARENT_SCOPE)
        list(APPEND UTOPIA_BUILD_INCLUDES ${EIGEN3_INCLUDE_DIR})
        # target_include_directories(utopia PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()



# TODO FIGURE THIS OUT.
# if(TARGET utopia_opencl)
#     target_link_libraries(utopia PUBLIC utopia_opencl)
# endif()

# if(TARGET utopia_petsc)
#     target_link_libraries(utopia PUBLIC utopia_petsc)
# endif()

# if(TARGET utopia_blas)
#     target_link_libraries(utopia PUBLIC utopia_blas)
# endif()

# if(TARGET utopia_cuda_cxx)
#     target_link_libraries(utopia PUBLIC utopia_cuda_cxx)
# endif()

# if(TARGET utopia_trilinos)
#     target_link_libraries(utopia PUBLIC utopia_trilinos)
# endif()

# if(TARGET utopia_m3elinsol)
#     target_link_libraries(utopia PUBLIC utopia_m3elinsol)
# endif()
