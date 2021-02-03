# InstallPetsc.cmake

if(NOT PETSC_FOUND)
    # git clone -b maint https://gitlab.com/petsc/petsc.git petsc
    include(ExternalProject)

    if(UTOPIA_DEPENDENCIES_DIR)
        set(PETSC_INSTALL_DIR ${UTOPIA_DEPENDENCIES_DIR}/petsc)
    else()
        set(PETSC_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../external/petsc)
    endif()

    set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
    set(PETSC_URL https://gitlab.com/petsc/petsc.git)
    set(PETSC_SOURCE_DIR ${STAGE_DIR}/petsc)
    set(PETSC_BIN_DIR ${STAGE_DIR}/petsc/bin)

    set(PETSC_MPI_BASE_DIR $ENV{MPI_DIR})
    set(MAKE_COMMAND "make")

    set(PETSC_CONFIG_ARGS $ENV{PETSC_CONFIG_ARGS})
    set(PETSC_CONFIG_ARGS
        ${PETSC_CONFIG_ARGS}
        --with-mpi=1
        --download-scalapack=yes
        # --with-scalapack=1 --with-scalapack-dir=/opt/local
        --download-hypre=yes
        --with-cxx-dialect=C++11
        -download-superlu_dist=yes
        --download-superlu=yes
        --download-mumps=yes
        -with-debugging=0)

    # ##########################################################################

    if(UTOPIA_ENABLE_PETSC_DM_PLEX)
        # DMPlex dependencies
        set(PETSC_CONFIG_ARGS
            ${PETSC_CONFIG_ARGS}
            --download-netcdf
            --download-pnetcdf
            --download-exodusii
            --download-zlib
            --download-triangle
            --download-ctetgen)

        if(HDF5_DIR)
            set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS}
                                  --with-hdf5-dir=${HDF5_DIR})
        else()
            set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS} --download-hdf5)
        endif()

        if(UTOPIA_ENABLE_CGNS)
            set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS} --with-cgns=yes
                                  --with-cgns-dir=/opt/local)
            message(STATUS "USING CGNS")
        endif()
    endif()

    # ##########################################################################

    ExternalProject_Add(
        petsc
        UPDATE_COMMAND "" # FIXME
        BUILD_IN_SOURCE 1
        PREFIX ${STAGE_DIR}
        GIT_REPOSITORY ${PETSC_URL}
        GIT_TAG maint
        DOWNLOAD_DIR ${STAGE_DIR}
        INSTALL_DIR ${PETSC_INSTALL_DIR}
        LOG_CONFIGURE 1
        LOG_BUILD 1
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
                          ${PETSC_CONFIG_ARGS}
        BUILD_COMMAND ${MAKE_COMMAND}
        INSTALL_COMMAND make install
        # COMMAND       ${MAKE_COMMAND}
    )

    set_target_properties(petsc PROPERTIES EXCLUDE_FROM_ALL TRUE)

    # set(PETSC_DIR ${PETSC_INSTALL_DIR})
endif()
