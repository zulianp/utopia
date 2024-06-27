# InstallTrilinos.cmake

if(NOT Trilinos_FOUND)
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

  if(UTOPIA_ENABLE_ENV_READ)
    set(HDF5_DIR $ENV{HDF5_DIR})
  endif()

  if(UTOPIA_ENABLE_CLUSTER)
    message(
      STATUS
        "On Cray System: Adding extra variables to find Netcdf, Pnetcdf and local install of SuperLU."
    )

    set(MPI_DIR $ENV{CRAY_MPICH_BASEDIR})

    list(
      APPEND
      TRILINOS_CMAKE_ARGS
      "-DNetcdf_INCLUDE_DIRS=$ENV{NETCDF_DIR}/include/;$ENV{PNETCDF_DIR}/include"
      "-DNetcdf_LIBRARY_DIRS=$ENV{NETCDF_DIR}/lib/;$ENV{PNETCDF_DIR}/lib"
      "-TPL_Netcdf_INCLUDE_DIRS=$ENV{NETCDF_DIR}/lib/;$ENV{PNETCDF_DIR}/lib"
      "-DSuperLU_INCLUDE_DIRS=$ENV{SuperLU_DIR}/include"
      "-DSuperLU_LIBRARY_DIRS=$ENV{SuperLU_DIR}/lib64")
  endif()

  list(
    APPEND
    TRILINOS_CMAKE_ARGS
    "-DCMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR}"
    "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    "-DCMAKE_CXX_STANDARD=17"
    "-DBUILD_SHARED_LIBS=OFF"
    "-DAmesos2_ENABLE_EXPLICIT_INSTANTIATION=ON"
    "-DIfpack2_ENABLE_EXPLICIT_INSTANTIATION=ON"
    "-DMueLu_ENABLE_EXPLICIT_INSTANTIATION=ON"
    "-DTPL_ENABLE_MPI=ON"
    # "-DMPI_BASE_DIR=${MPI_DIR}"
    "-DTPL_ENABLE_Netcdf:BOOL=ON"
    "-DTPL_ENABLE_Pnetcdf:BOOL=OFF"
    "-DTeuchos_ENABLE_EXPLICIT_INSTANTIATION=ON"
    "-DTpetraCore_ENABLE_EXPLICIT_INSTANTIATION=ON"
    "-DTpetra_ENABLE_EXPLICIT_INSTANTIATION=ON"
    "-DTpetra_INST_COMPLEX_DOUBLE=OFF"
    "-DTpetra_INST_DOUBLE:BOOL=ON"
    "-DTpetra_INST_INT_LONG:BOOL=ON"
    "-DTpetra_INST_INT_LONG_LONG:BOOL=OFF"
    "-DTrilinos_ALLOW_NO_PACKAGES:BOOL=OFF"
    "-DTrilinos_ASSERT_DEFINED_DEPENDENCIES=FATAL_ERROR"
    "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF "
    "-DTrilinos_ENABLE_Amesos2:BOOL=ON"
    "-DTrilinos_ENABLE_AztecOO:BOOL=OFF"
    "-DTrilinos_ENABLE_Belos:BOOL=ON"
    "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
    "-DTrilinos_ENABLE_Epetra:BOOL=OFF"
    "-DTrilinos_ENABLE_EpetraExt:BOOL=OFF"
    "-DTrilinos_ENABLE_Gtest:BOOL=OFF"
    "-DTrilinos_ENABLE_Ifpack2:BOOL=ON"
    "-DTrilinos_ENABLE_Intrepid2:BOOL=ON"
    "-DTrilinos_ENABLE_MueLu:BOOL=ON"
    "-DTrilinos_ENABLE_NOX=ON "
    "-DTrilinos_ENABLE_SEACASEpu:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASExodiff:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASExodus:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASIoss:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASNemslice:BOOL=ON"
    "-DTrilinos_ENABLE_SEACASNemspread:BOOL=ON"
    "-DTrilinos_ENABLE_STKBalance:BOOL=OFF"
    "-DTrilinos_ENABLE_STKIO:BOOL=ON"
    "-DTrilinos_ENABLE_STKMesh:BOOL=ON"
    "-DTrilinos_ENABLE_STKSearch:BOOL=ON"
    "-DTrilinos_ENABLE_STKSimd:BOOL=ON"
    "-DTrilinos_ENABLE_STKTopology:BOOL=ON"
    "-DTrilinos_ENABLE_STKTransfer:BOOL=ON"
    "-DTrilinos_ENABLE_STKUnit_test_utils:BOOL=OFF"
    "-DTrilinos_ENABLE_STKUnit_tests:BOOL=OFF"
    "-DTrilinos_ENABLE_STKUtil:BOOL=ON"
    "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
    "-DTrilinos_ENABLE_Tpetra:BOOL=ON"
    "-DTrilinos_ENABLE_TpetraCore=ON"
    "-DTrilinos_ENABLE_Zoltan2:BOOL=ON"
    "-DTrilinos_ENABLE_Zoltan:BOOL=ON"
    "-DTPL_ENABLE_SuperLU:BOOL=OFF"
    # "-DTPL_ENABLE_SuperLU:BOOL=ON"
    "-DXpetra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
    "-DTrilinos_ENABLE_Percept:BOOL=ON"
    "-DTPL_ENABLE_HDF5:BOOL=ON"
    "-DXpetra_ENABLE_EXPLICIT_INSTANTIATION=ON"
    "-DTrilinos_ENABLE_Kokkos=ON"
    "-DTPL_ENABLE_HDF5=ON"
    "-DHDF5_INCLUDE_DIRS=${HDF5_DIR}/include/"
    "-DHDF5_LIBRARY_DIRS=${HDF5_DIR}/lib/"
    "-DTrilinos_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR=${CMAKE_SOURCE_DIR}/../external/"
    "-DTrilinos_ENABLE_EXAMPLES=OFF")

  # For cuda
  if(UTOPIA_ENABLE_GPU)
    list(
      APPEND
      TRILINOS_CMAKE_ARGS
      "-DCMAKE_CXX_COMPILER=${STAGE_DIR}/src/trilinos/packages/kokkos/bin/nvcc_wrapper"
      "-DKokkos_ENABLE_CUDA=ON"
      "-DKokkos_ENABLE_CUDA_CONSTEXPR=ON"
      "-DKokkos_ENABLE_CUDA_LAMBDA=ON"
      "-DCMAKE_CXX_STANDARD=17"
      "-DKokkos_ARCH_PASCAL61=ON "
      "-DTpetra_INST_CUDA=ON")
  endif()

  if(EXISTS ${HDF5_DIR})
    ExternalProject_Add(
      trilinos
      UPDATE_COMMAND "" # FIXME
      PREFIX ${STAGE_DIR}
      GIT_REPOSITORY ${TRILINOS_URL}
      # GIT_TAG trilinos-release-15-0-0
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
    set(Trilinos_DIR ${TRILINOS_INSTALL_DIR})
  else()
    message(
      FATAL_ERROR
        "Please set the following variables for trilinos to install correctly:\nHDF5_DIR: Folder where HDF is located."
    )
  endif()
endif()
