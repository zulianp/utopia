cmake .. \
-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/Trilinos \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_STANDARD=14 \
-DBUILD_SHARED_LIBS=OFF \
-DAmesos2_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DIfpack2_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DMueLu_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTPL_ENABLE_MPI=ON \
-DTPL_ENABLE_Netcdf:BOOL=ON \
-DTPL_ENABLE_Pnetcdf=OFF \
-DTeuchos_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTpetraCore_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTpetra_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTpetra_INST_COMPLEX_DOUBLE=OFF \
-DTpetra_INST_DOUBLE:BOOL=ON \
-DTpetra_INST_INT_LONG:BOOL=ON \
-DTpetra_INST_INT_LONG_LONG:BOOL=OFF \
-DTrilinos_ALLOW_NO_PACKAGES:BOOL=OFF \
-DTrilinos_ASSERT_MISSING_PACKAGES=ON \
-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF  \
-DTrilinos_ENABLE_Amesos2:BOOL=ON \
-DTrilinos_ENABLE_AztecOO:BOOL=OFF \
-DTrilinos_ENABLE_Belos:BOOL=ON \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTrilinos_ENABLE_Epetra:BOOL=OFF \
-DTrilinos_ENABLE_EpetraExt:BOOL=OFF \
-DTrilinos_ENABLE_Gtest:BOOL=OFF \
-DTrilinos_ENABLE_Ifpack2:BOOL=ON \
-DTrilinos_ENABLE_Intrepid2:BOOL=ON \
-DTrilinos_ENABLE_MueLu:BOOL=ON \
-DTrilinos_ENABLE_NOX=ON  \
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
-DTrilinos_ENABLE_STKUnit_test_utils:BOOL=OFF \
-DTrilinos_ENABLE_STKUnit_tests:BOOL=OFF \
-DTrilinos_ENABLE_STKUtil:BOOL=ON \
-DTrilinos_ENABLE_TESTS:BOOL=OFF \
-DTrilinos_ENABLE_Tpetra:BOOL=ON \
-DTrilinos_ENABLE_TpetraCore=ON \
-DTrilinos_ENABLE_Zoltan2:BOOL=ON \
-DTrilinos_ENABLE_Zoltan:BOOL=ON \
-DTPL_ENABLE_SuperLU:BOOL=OFF \
-DXpetra_ENABLE_EXPLICIT_INSTANTIATION=ON
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
