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
        message(STATUS "TRILINOS_MPI_BASE_DIR=${TRILINOS_MPI_BASE_DIR}")
        message(STATUS "MPI_CXX_COMPILER=${MPI_CXX_COMPILER}")
        message(STATUS "MPI_C_COMPILER=${MPI_C_COMPILER}")
    endif()

    set(TRILINOS_CXX_COMPILER ${CMAKE_CXX_COMPILER})
    set(TRILINOS_C_COMPILER ${CMAKE_C_COMPILER})

    if(MPI_CXX_COMPILER)
        list(APPEND TRILINOS_CMAKE_ARGS "-DMPI_USE_COMPILER_WRAPPERS=ON"
             "-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}")
    endif()

    if(MPI_C_COMPILER)
        list(APPEND TRILINOS_CMAKE_ARGS "-DMPI_C_COMPILER=${MPI_C_COMPILER}")
    endif()

    if(MPI_Fortan_COMPILER)
        list(APPEND TRILINOS_CMAKE_ARGS
             "-DMPI_Fortan_COMPILER=${MPI_Fortan_COMPILER}")
    endif()

    list(
        APPEND
        TRILINOS_CMAKE_ARGS
        "-DCMAKE_CXX_STANDARD=17"
        "-DCMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR}"
        "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        "-DCMAKE_CXX_COMPILER=${TRILINOS_CXX_COMPILER}"
        "-DCMAKE_C_COMPILER=${TRILINOS_C_COMPILER}"
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
        "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
        "-DTpetra_INST_DOUBLE:BOOL=ON"
        "-DTpetra_INST_INT_LONG:BOOL=ON"
        "-DTpetra_INST_INT_LONG_LONG:BOOL=OFF"
        "-DTpetra_INST_COMPLEX_DOUBLE=OFF"
        "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
        "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF "
        "-DAmesos2_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DIfpack2_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DMueLu_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DRTOp_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DStratimikos_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTeuchos_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DThyra_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTpetraCore_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTpetra_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON"
        "-DXpetra_ENABLE_EXPLICIT_INSTANTIATION=ON")

    list(
        APPEND
        TRILINOS_CMAKE_ARGS
        "-DTrilinos_ENABLE_STKMesh:BOOL=ON"
        "-DTrilinos_ENABLE_STKSimd:BOOL=ON"
        "-DTrilinos_ENABLE_STKIO:BOOL=ON"
        "-DTrilinos_ENABLE_STKTransfer:BOOL=ON"
        "-DTrilinos_ENABLE_STKSearch:BOOL=ON"
        "-DTrilinos_ENABLE_STKUtil:BOOL=ON"
        "-DTrilinos_ENABLE_STKTopology:BOOL=ON"
        "-DTrilinos_ENABLE_STKBalance:BOOL=OFF"
        "-DTrilinos_ENABLE_STKUnit_tests:BOOL=OFF"
        "-DTrilinos_ENABLE_STKUnit_test_utils:BOOL=OFF"
        "-DTrilinos_ENABLE_Gtest:BOOL=ON"
        "-DTrilinos_ENABLE_SEACASExodus:BOOL=ON"
        "-DTrilinos_ENABLE_SEACASEpu:BOOL=ON"
        "-DTrilinos_ENABLE_SEACASExodiff:BOOL=ON"
        "-DTrilinos_ENABLE_SEACASNemspread:BOOL=ON"
        "-DTrilinos_ENABLE_SEACASNemslice:BOOL=ON"
        "-DTrilinos_ENABLE_SEACASIoss:BOOL=ON"
        "-DTrilinos_ENABLE_Percept:BOOL=ON")

    list(APPEND TRILINOS_CMAKE_ARGS "-DTrilinos_ENABLE_Intrepid2:BOOL=ON")

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
