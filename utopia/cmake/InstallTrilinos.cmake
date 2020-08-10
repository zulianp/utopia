# InstallTrilinos.cmake

if(NOT TRILINOS_FOUND)
    # git clone https://github.com/trilinos/Trilinos.git
    include(ExternalProject)

    if(UTOPIA_DEPENDENCIES_DIR)
        set(TRILINOS_INSTALL_DIR ${UTOPIA_DEPENDENCIES_DIR}/Trilinos)
    else()
        set(TRILINOS_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../external/Trilinos)
    endif()

    set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
    set(TRILINOS_URL https://github.com/trilinos/Trilinos.git)
    set(TRILINOS_SOURCE_DIR ${STAGE_DIR}/Trilinos)
    set(TRILINOS_BIN_DIR ${STAGE_DIR}/Trilinos/bin)

    if(NOT MPI_DIR)
        set(MPI_DIR $ENV{MPI_DIR})
    endif()

    if(MPI_DIR)
        set(TRILINOS_MPI_BASE_DIR ${MPI_DIR})
    endif()

    list(
        APPEND
        TRILINOS_CMAKE_ARGS
        "-DCMAKE_CXX_STANDARD=14"
        "-DCMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR}"
        "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        "-DENABLE_SANITIZER=${ENABLE_SANITIZER}"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
        "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        "-DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}"
        "-DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}"
        "-DTPL_ENABLE_MPI=ON"
        "-DMPI_BASE_DIR=${TRILINOS_MPI_BASE_DIR}"
        "-DTrilinos_ENABLE_Tpetra=ON"
        "-DTrilinos_ENABLE_TpetraCore=ON"
        "-DTrilinos_ENABLE_Belos=ON"
        "-DTrilinos_ENABLE_Amesos2=ON"
        "-DTrilinos_ENABLE_Ifpack2=ON"
        "-DTrilinos_ENABLE_MueLu=ON"
        "-DTrilinos_ENABLE_NOX=ON "
        "-DAmesos2_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DIfpack2_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DKOKKOS_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DMueLu_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DRTOp_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DStratimikos_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTeuchos_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DThyra_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTpetraCore_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTpetra_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DXpetra_ENABLE_EXPLICIT_INSTANTIATION=ON")

    # For cuda
    if(UTOPIA_ENABLE_CUDA)
        list(
            APPEND
            TRILINOS_CMAKE_ARGS
            "-DCMAKE_CXX_COMPILER=/home/zulian/bin/bin/nvcc_wrapper"
            "-DKokkos_ENABLE_CUDA=ON"
            "-DKokkos_ENABLE_CUDA_CONSTEXPR=ON"
            "-DKokkos_ENABLE_CUDA_LAMBDA=ON"
            "-DKokkos_ENABLE_CUDA_UVM=ON"
            "-DCMAKE_CXX_STANDARD=14"
            "-DKokkos_ARCH_PASCAL61=ON "
            "-DTpetra_INST_CUDA=ON")
    endif()

    ExternalProject_Add(
        trilinos
        UPDATE_COMMAND "" # FIXME
        PREFIX ${STAGE_DIR}
        GIT_REPOSITORY ${TRILINOS_URL}
        DOWNLOAD_DIR ${STAGE_DIR}
        INSTALL_DIR ${TRILINOS_INSTALL_DIR}
        # BINARY_DIR                      ${TRILINOS_SOURCE_DIR}
        CMAKE_ARGS "${TRILINOS_CMAKE_ARGS}"
        LOG_CONFIGURE 1
        LOG_BUILD 1
        BUILD_COMMAND ${CMAKE_COMMAND} -E echo "Starting $<CONFIG> build"
        COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config $<CONFIG>
        COMMAND ${CMAKE_COMMAND} -E echo "$<CONFIG> build complete")

    set_target_properties(trilinos PROPERTIES EXCLUDE_FROM_ALL TRUE)

    # set(TRILINOS_DIR ${TRILINOS_INSTALL_DIR}) set(TRILINOS_DIR
    # ${TRILINOS_INSTALL_DIR} PARENT_SCOPE) set(TRILINOS_FOUND TRUE)

endif()
