utopia_add_library(utopia_interoperability ${CMAKE_CURRENT_SOURCE_DIR} "." "utopia_core;utopia_utils;utopia_solvers")

# if(TARGET utopia_opencl)
# 	target_link_libraries(utopia_interoperability PUBLIC utopia_opencl)
# endif()

if(TARGET utopia_petsc)
    target_link_libraries(utopia_interoperability PUBLIC utopia_petsc)
endif()

# if(TARGET utopia_blas)
# 	target_link_libraries(utopia_interoperability PUBLIC utopia_blas)
# endif()

# if(TARGET utopia_cuda_cxx)
# 	target_link_libraries(utopia_interoperability PUBLIC utopia_cuda_cxx)
# endif()

if(TARGET utopia_trilinos)
    target_link_libraries(utopia_interoperability PUBLIC utopia_trilinos)
endif()

# if(TARGET utopia_m3elinsol)
# 	target_link_libraries(utopia_interoperability utopia_m3elinsol)
# endif()

