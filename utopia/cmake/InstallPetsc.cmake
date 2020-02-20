# InstallTrilinos.cmake

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

    set(PETSC_CONFIG_ARGS "${PETSC_CONFIG_ARGS}")

    ExternalProject_Add(
        petsc
        UPDATE_COMMAND      "" #FIXME
        BUILD_IN_SOURCE     1
        PREFIX              ${STAGE_DIR}
        GIT_REPOSITORY      ${PETSC_URL}
        GIT_TAG             v3.11.2
        DOWNLOAD_DIR        ${STAGE_DIR}
        INSTALL_DIR         ${PETSC_INSTALL_DIR}
        LOG_CONFIGURE       1
        LOG_BUILD           1
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --with-mpi=1 --download-superlu_dist=yes --download-superlu=yes --download-mumps=yes --download-scalapack=yes --download-hypre=yes --with-cxx-dialect=C++11 --with-debugging=0 #-with-mpi-include=${MPI_INCLUDE_PATH} --with-mpi-lib=${MPI_LIBRARIES}
        BUILD_COMMAND ${MAKE_COMMAND}
        INSTALL_COMMAND     make install
        # COMMAND       ${MAKE_COMMAND}
    )

    set(PETSC_DIR ${PETSC_INSTALL_DIR})
    set(PETSC_DIR ${PETSC_INSTALL_DIR} PARENT_SCOPE)
endif()
