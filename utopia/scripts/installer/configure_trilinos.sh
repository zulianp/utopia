# Tested on Piz Daint (With master, fmt folder needs to be copied to the installation folder Trilinos/include)
# https://github.com/trilinos/Trilinos.git
# commit 77005adad6d625dbf62009620ffdc4ffa06b9fac

cmake .. \
-DCMAKE_INSTALL_PREFIX=$SCRATCH/installations/Trilinos \
-DNetcdf_LIBRARY_DIRS=$NETCDF_DIR/lib/   \
-DTPL_Netcdf_INCLUDE_DIRS=$NETCDF_DIR/include/ \
-DTPL_Netcdf_LIBRARIES=$NETCDF_DIR/lib/libnetcdf.a \
-DCMAKE_C_COMPILER=cc \
-DCMAKE_CXX_COMPILER=CC \
-DTrilinos_ENABLE_Fortran:BOOL=OFF \
-DAmesos2_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DBUILD_SHARED_LIBS=OFF \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_STANDARD=17 \
-DTeuchos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTpetra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTpetra_INST_COMPLEX_DOUBLE=OFF \
-DTpetra_INST_DOUBLE:BOOL=ON \
-DTpetra_INST_INT_LONG:BOOL=ON \
-DTpetra_INST_INT_LONG_LONG:BOOL=OFF \
-DTpetraCore_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTPL_ENABLE_Boost:BOOL=OFF \
-DTPL_ENABLE_HDF5:BOOL=OFF \
-DTPL_ENABLE_MPI:BOOL=ON \
-DTPL_ENABLE_Netcdf:BOOL=ON \
-DTPL_ENABLE_Pnetcdf:BOOL=OFF \
-DTPL_ENABLE_SuperLU:BOOL=OFF \
-DTrilinos_ALLOW_NO_PACKAGES:BOOL=OFF \
-DTrilinos_ASSERT_DEFINED_DEPENDENCIES=FATAL_ERROR \
-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF  \
-DTrilinos_ENABLE_Amesos2:BOOL=ON \
-DTrilinos_ENABLE_AztecOO:BOOL=OFF \
-DTrilinos_ENABLE_Belos:BOOL=ON \
-DTrilinos_ENABLE_Epetra:BOOL=OFF \
-DTrilinos_ENABLE_EpetraExt:BOOL=OFF \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTrilinos_ENABLE_Gtest:BOOL=OFF \
-DTrilinos_ENABLE_Ifpack2:BOOL=ON \
-DIfpack2_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTrilinos_ENABLE_Intrepid2:BOOL=ON \
-DTrilinos_ENABLE_Kokkos=ON  \
-DTrilinos_ENABLE_MueLu:BOOL=OFF \
-DMueLu_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
-DTrilinos_ENABLE_NOX:BOOL=OFF  \
-DTrilinos_ENABLE_Percept:BOOL=ON \
-DTrilinos_ENABLE_SEACASEpu:BOOL=ON \
-DTrilinos_ENABLE_SEACASExodiff:BOOL=ON \
-DTrilinos_ENABLE_SEACASExodus:BOOL=ON \
-DTrilinos_ENABLE_SEACASIoss:BOOL=ON \
-DTrilinos_ENABLE_SEACASNemslice:BOOL=ON \
-DTrilinos_ENABLE_SEACASNemspread:BOOL=ON \
-DTrilinos_ENABLE_STKBalance:BOOL=OFF \
-DTrilinos_ENABLE_STKIO:BOOL=ON \
-DTrilinos_ENABLE_STKMesh:BOOL=ON \
-DTrilinos_ENABLE_STKSearch:BOOL=ON \
-DTrilinos_ENABLE_STKSimd:BOOL=ON \
-DTrilinos_ENABLE_STKTopology:BOOL=ON \
-DTrilinos_ENABLE_STKTransfer:BOOL=ON \
-DTrilinos_ENABLE_STKUnit_tests:BOOL=OFF \
-DTrilinos_ENABLE_STKUtil:BOOL=ON \
-DTrilinos_ENABLE_TESTS:BOOL=OFF \
-DTrilinos_ENABLE_STKUnit_test_utils:BOOL=OFF \
-DTrilinos_ENABLE_Tpetra:BOOL=ON \
-DTrilinos_ENABLE_TpetraCore:BOOL=ON \
-DTrilinos_ENABLE_Zoltan2:BOOL=ON \
-DTrilinos_ENABLE_Zoltan:BOOL=ON \
-DXpetra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON  && \
make -j40 && make install && \
cp -r ../packages/seacas/libraries/ioss/src/private_copy_fmt/fmt/ $SCRATCH/installations/Trilinos/include/fmt


# Enable and specify HDF5 if linkage errors.  
# Trilinos compiled no problems afterwards.

#  Warning, Trilinos_ASSERT_MISSING_PACKAGES='ON' is set and is no longer
#  supported! Please set Trilinos_ASSERT_DEFINED_DEPENDENCIES instead

# -DTPL_ENABLE_SuperLU=ON \
# -DTPL_ENABLE_Boost:BOOL=ON \
# -DTPL_ENABLE_HDF5=ON \
# -DTPL_HDF5_INCLUDE_DIRS=$HDF5_INC_DIR \
# -DTPL_HDF5_LIBRARY_DIRS=$HDF5_LIB_DIR \
# -DTPL_Netcdf_INCLUDE_DIRS=$NETCDF_INC_DIR \
# -DTPL_SuperLU_INCLUDE_DIRS=$SUPERLU_INC_DIR \
# -DBoost_INCLUDE_DIRS:PATH=$BOOST_INC_DIR \
# -DMPI_BASE_DIR:PATH=$MPI_BASE_DIR \
# -DNetcdf_LIBRARIES=$NETCDF_LIB_DIR \
