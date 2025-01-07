# InstallPetsc.cmake
# If it does not work try directly in the terminal
# git clone -b main https://gitlab.com/petsc/petsc.git petsc
# ./configure --prefix=$INSTALL_DIR/petsc --with-mpi=1 --download-scalapack=yes --download-hypre=yes --download-metis=yes --download-parmetis=yes --download-mumps=yes --with-debugging=0

if(NOT CYGWIN)
  if(NOT PETSC_FOUND)
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
        --download-hypre=yes
        --download-metis=yes
        --download-parmetis=yes
        # --with-cxx-dialect=C++11
        --download-mumps=yes
        --with-debugging=0)

    if(UTOPIA_ENABLE_SLEPC)
      set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS} --download-slepc=yes)
    endif()

    if(UTOPIA_PETSC_ENABLE_SUPERLU)
      set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS} --download-superlu_dist=yes
                            --download-superlu=yes)
    endif()

    if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_HOST_SYSTEM_VERSION
                                                    VERSION_GREATER_EQUAL 20)
      set(APPLE_GREATER_OR_EQUAL_THAN_BIG_SUR TRUE)
    endif()

    if(NOT APPLE_GREATER_OR_EQUAL_THAN_BIG_SUR AND NOT UNIX)
      if(BLAS_LIBRARIES)
        message(STATUS "[InstallPetsc.cmake] BLAS_LIBRARIES=${BLAS_LIBRARIES}")
        set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS}
                              --with-blas-lib=${BLAS_LIBRARIES})
      endif()

      if(LAPACK_LIBRARIES)
        message(
          STATUS "[InstallPetsc.cmake] LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
        set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS}
                              --with-lapack-lib=${LAPACK_LIBRARIES})
      endif()
    elseif(CrayLinuxEnvironment)
      message(STATUS "Running on Cray System, Skipping blas and lapack pointers for PetscInstall.cmake")
    elseif(Linux)
      message(STATUS "Running on Linux System, Skipping blas and lapack pointers for PetscInstall.cmake")
    else()
      message(
        STATUS
          "Running on MacOS >= Big Sur, Skipping blas and lapack pointers for PetscInstall.cmake"
      )
    endif()
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
        set(PETSC_CONFIG_ARGS ${PETSC_CONFIG_ARGS} --with-hdf5-dir=${HDF5_DIR})
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

    if(DEFINED $ENV{PETSC_DIR})
      message(FATAL_ERROR "Please unset PETSC_DIR, PETSC_ARCH to enable local petsc target install.")
    else()
      ExternalProject_Add(
        petsc
        UPDATE_COMMAND "" # FIXME
        BUILD_IN_SOURCE 1
        PREFIX ${STAGE_DIR}
        GIT_REPOSITORY ${PETSC_URL}
        GIT_TAG v3.20.2
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
      set($ENV{PETSC_DIR} "${PETSC_INSTALL_DIR}")

    endif()


    # ${PETSC_INSTALL_DIR})

  endif()
endif()
