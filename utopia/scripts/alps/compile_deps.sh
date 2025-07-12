#!/usr/bin/env bash

set -e
set -x

export DEPS=$HOME/deps
export INSTALL_DIR=$DEPS/installations

mkdir -p $DEPS
cd $DEPS

####################
# Install moonolith
####################

git clone https://bitbucket.org/zulianp/par_moonolith.git
cd $DEPS/par_moonolith

mkdir -p build && cd build
cmake .. \  
		 -DCMAKE_INSTALL_PREFIX=$DEPS/installations/par_moonolith \
		 -DCMAKE_CXX_FLAGS='-mcpu=neoverse-v2 -O3 -ffast-math'

make -j72

####################
# Install petsc
####################

cd $DEPS
git clone https://github.com/petsc/petsc.git
cd $DEPS/petsc

./configure CFLAGS='-mcpu=neoverse-v2 -O3 -ffast-math -funroll-loops -march=armv9-a' CXXFLAGS='-mcpu=neoverse-v2 -O3 -ffast-math' --prefix=$DEPS/installations/petsc --with-mpi=1 --download-scalapack=yes --download-hypre=yes --download-metis=yes --download-parmetis=yes --with-cxx-dialect=C++11 --download-mumps=yes --with-debugging=0 --download-superlu_dist=yes --download-superlu=yes  --download-netcdf --download-pnetcdf --download-exodusii --download-zlib --download-triangle --download-ctetgen --download-hdf5  --download-fblaslapack=1 

make PETSC_DIR=/users/zulianp/deps/petsc PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=/users/zulianp/deps/petsc PETSC_ARCH=arch-linux-c-opt install

####################
# Install Trilinos
####################

cd $DEPS
git clone https://github.com/trilinos/Trilinos.git
cd $DEPS/Trilinos

mkdir -p build && cd build

cmake .. \
-DCMAKE_INSTALL_PREFIX=$DEPS/installations/Trilinos \
-DNetcdf_LIBRARY_DIRS=$DEPS/installations/petsc/lib/   \
-DTPL_Netcdf_INCLUDE_DIRS=$DEPS/installations/petsc/include/ \
-DTPL_Netcdf_LIBRARIES=$DEPS/installations/petsc/lib/libnetcdf.a \
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
-DXpetra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON

make -j72 && make install
 # && \
# cp -r ../packages/seacas/libraries/ioss/src/private_copy_fmt/fmt/ $SCRATCH/installations/Trilinos/include/fmt

####################
# Install Utopia
####################

cd $HOME
git clone https://github.com/zulianp/utopia.git
git submodule update --init --recursive
cd $HOME/utopia/utopia
mkdir build && cd build

cmake .. \
	-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia \
	-DUTOPIA_ENABLE_FLUYA_MODE=ON \
	-DUTOPIA_INSTALL_YAML=ON \
	-DPetsc_DIR=$INSTALL_DIR/petsc \ 
	-DTrilinos_DIR=$INSTALL_DIR/Trilinos/lib/cmake/Trilinos 

make yaml-cpp

cmake ..

make -j72 && make install
