# # InstallPetsc.cmake

# if(NOT PETSC_FOUND)

#     #  git clone -b maint https://gitlab.com/petsc/petsc.git petsc
#     include(ExternalProject)

#     set(STAGE_DIR            "${CMAKE_BINARY_DIR}/stage")
#     set(PETSC_URL            https://gitlab.com/petsc/petsc.git)
#     set(PETSC_SOURCE_DIR     ${STAGE_DIR}/petsc)
#     set(PETSC_BIN_DIR        ${STAGE_DIR}/petsc/bin)
#     set(PETSC_INSTALL_DIR    ${CMAKE_SOURCE_DIR}/external/petsc)
#     set(PETSC_MPI_BASE_DIR   $ENV{MPI_DIR})
#     set(MAKE_COMMAND "make")

#     set(PETSC_CONFIG_ARGS "$ENV{PETSC_CONFIG_ARGS}")

#     ExternalProject_Add(
#         petsc
#         UPDATE_COMMAND      "" #FIXME
#         BUILD_IN_SOURCE     1
#         PREFIX              ${STAGE_DIR}
#         GIT_REPOSITORY      ${PETSC_URL}
#         #GIT_TAG             tags/v3.11.2
#         GIT_TAG		        maint
# 	    DOWNLOAD_DIR        ${STAGE_DIR}
#         INSTALL_DIR         ${PETSC_INSTALL_DIR}
#         LOG_CONFIGURE       1
#         LOG_BUILD           1
#         CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --with-mpi=1 --download-scalapack=yes --download-hypre=yes --with-cxx-dialect=C++11 --with-debugging=0  --with-64-bit-indices=1 ${PETSC_CONFIG_ARGS} #-with-mpi-include=${MPI_INCLUDE_PATH} --with-mpi-lib=${MPI_LIBRARIES} #--download-superlu_dist=yes --download-superlu=yes --download-mumps=yes
#         BUILD_COMMAND       ${MAKE_COMMAND}
#         INSTALL_COMMAND     make install
#         # COMMAND       ${MAKE_COMMAND}
#     )

#     set_target_properties(petsc PROPERTIES EXCLUDE_FROM_ALL TRUE)

#     # set(PETSC_DIR ${PETSC_INSTALL_DIR})
# endif()


# InstallPetsc.cmake

if(NOT PETSC_FOUND)

    #  git clone -b maint https://gitlab.com/petsc/petsc.git petsc
    include(ExternalProject)

    set(STAGE_DIR            "${CMAKE_BINARY_DIR}/stage")
    set(PETSC_URL            https://gitlab.com/petsc/petsc.git)
    set(PETSC_SOURCE_DIR     ${STAGE_DIR}/petsc)
    set(PETSC_BIN_DIR        ${STAGE_DIR}/petsc/bin)
    set(PETSC_INSTALL_DIR    ${CMAKE_SOURCE_DIR}/external/petsc)
    set(PETSC_MPI_BASE_DIR   $ENV{MPI_DIR})
    set(MAKE_COMMAND "make")

    # set(PETSC_CONFIG_ARGS "$ENV{PETSC_CONFIG_ARGS} ${PETSC_CONFIG_ARGS}")
    set(PETSC_CONFIG_ARGS "$ENV{PETSC_CONFIG_ARGS}")

    if(UTOPIA_ENABLE_CGNS)
        set(PETSC_CONFIG_ARGS "--with-cgns=yes --with-cgns-dir=/opt/local")
        message(STATUS "USING CGNS")
    endif()

    ExternalProject_Add(
        petsc
        UPDATE_COMMAND      "" #FIXME
        BUILD_IN_SOURCE     1
        PREFIX              ${STAGE_DIR}
        GIT_REPOSITORY      ${PETSC_URL}
        GIT_TAG             maint
        DOWNLOAD_DIR        ${STAGE_DIR}
        INSTALL_DIR         ${PETSC_INSTALL_DIR}
        LOG_CONFIGURE       1
        LOG_BUILD           1
        CONFIGURE_COMMAND   <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --with-mpi=1 --download-scalapack=yes --download-hypre=yes --with-cxx-dialect=C++11 -download-superlu_dist=yes --download-superlu=yes --download-mumps=yes --download-hdf5 --download-netcdf --download-pnetcdf --download-exodusii --download-zlib --download-triangle --download-ctetgen -with-debugging=0 ${PETSC_CONFIG_ARGS}  #-with-mpi-include=${MPI_INCLUDE_PATH} --with-mpi-lib=${MPI_LIBRARIES} -
        BUILD_COMMAND       ${MAKE_COMMAND}
        INSTALL_COMMAND     make install
        # COMMAND       ${MAKE_COMMAND}
    )

    set_target_properties(petsc PROPERTIES EXCLUDE_FROM_ALL TRUE)

    # set(PETSC_DIR ${PETSC_INSTALL_DIR})
endif()
