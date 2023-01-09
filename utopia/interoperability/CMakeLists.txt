scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "."
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


# utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")

# # if(TARGET utopia_opencl)
# # 	target_link_libraries(utopia_interoperability PUBLIC utopia_opencl)
# # endif()

# if(TARGET utopia_petsc)
#     target_link_libraries(utopia_interoperability PUBLIC utopia_petsc)
# endif()

# # if(TARGET utopia_blas)
# # 	target_link_libraries(utopia_interoperability PUBLIC utopia_blas)
# # endif()

# # if(TARGET utopia_cuda_cxx)
# # 	target_link_libraries(utopia_interoperability PUBLIC utopia_cuda_cxx)
# # endif()

# if(TARGET utopia_trilinos)
#     target_link_libraries(utopia_interoperability PUBLIC utopia_trilinos)
# endif()

# if(TARGET utopia_m3elinsol)
# 	target_link_libraries(utopia_interoperability utopia_m3elinsol)
# endif()

