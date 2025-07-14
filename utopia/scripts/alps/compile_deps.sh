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
# commit 0337bfe0b9dcc77abc5d44df0b7f57cdcdf2ff74
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
# FIXME!
git checkout c6f6bea
#ATTENTION!!! Manually edit Tpetra_SUPPORTED_KOKKOS_VERSION "4.3.1" to Tpetra_SUPPORTED_KOKKOS_VERSION "4.4.1" (check the version installed in the uenv)
# vim ./packages/tpetra/CMakeLists.txt
mkdir -p build && cd build

cmake .. \
-DCMAKE_INSTALL_PREFIX=$DEPS/installations/Trilinos \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_STANDARD=17 \
-DTPL_ENABLE_Kokkos=ON  -DTPL_Kokkos_ROOT=/user-environment/env/default \
-DTpetra_INST_DOUBLE:BOOL=ON \
-DTpetra_INST_INT_LONG:BOOL=ON \
-DTPL_ENABLE_MPI:BOOL=ON \
-DAmesos2_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DIfpack2_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTeuchos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTpetra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DXpetra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTpetraCore_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTrilinos_ASSERT_DEFINED_DEPENDENCIES=FATAL_ERROR \
-DTrilinos_ENABLE_Amesos2:BOOL=ON \
-DTrilinos_ENABLE_Belos:BOOL=ON \
-DTrilinos_ENABLE_Ifpack2:BOOL=ON \
-DTrilinos_ENABLE_Intrepid2:BOOL=ON \
-DTrilinos_ENABLE_Kokkos=ON  \
-DTrilinos_ENABLE_Percept:BOOL=ON \
-DTrilinos_ENABLE_SEACASEpu:BOOL=ON \
-DTrilinos_ENABLE_SEACASExodiff:BOOL=ON \
-DTrilinos_ENABLE_SEACASExodus:BOOL=ON \
-DTrilinos_ENABLE_SEACASIoss:BOOL=ON \
-DTrilinos_ENABLE_SEACASNemslice:BOOL=ON \
-DTrilinos_ENABLE_SEACASNemspread:BOOL=ON \
-DTrilinos_ENABLE_STKIO:BOOL=ON \
-DTrilinos_ENABLE_STKMesh:BOOL=ON \
-DTrilinos_ENABLE_STKSearch:BOOL=ON \
-DTrilinos_ENABLE_STKBalance:BOOL=ON \
-DTrilinos_ENABLE_STKSimd:BOOL=ON \
-DTrilinos_ENABLE_STKTopology:BOOL=ON \
-DTrilinos_ENABLE_STKTransfer:BOOL=ON \
-DTrilinos_ENABLE_STKUtil:BOOL=ON \
-DTrilinos_ENABLE_Tpetra:BOOL=ON \
-DTrilinos_ENABLE_TpetraCore:BOOL=ON \
-DTrilinos_ENABLE_Zoltan2:BOOL=ON \
-DTrilinos_ENABLE_Zoltan:BOOL=ON \
-DTPL_ENABLE_HDF5:BOOL=ON \
-DHDF5_INCLUDE_DIRS=$DEPS/installations/petsc/include/ \
-DHDF5_LIBRARY_DIRS=$DEPS/installations/petsc/lib/ \
-DTPL_HDF5_LIBRARIES=$DEPS/installations/petsc/lib/libhdf5.so \
-DTPL_ENABLE_Netcdf:BOOL=ON \
-DNetcdf_LIBRARY_DIRS=$DEPS/installations/petsc/lib/   \
-DTPL_Netcdf_INCLUDE_DIRS=$DEPS/installations/petsc/include/ \
-DTPL_ENABLE_Pnetcdf::BOOL=ON \
-DTPL_Netcdf_LIBRARIES="$DEPS/installations/petsc/lib/libpnetcdf.so;$DEPS/installations/petsc/lib/libnetcdf.so;-ldl" \
-DPnetcdf_INCLUDE_DIRS=$DEPS/installations/petsc/include \
-DBUILD_SHARED_LIBS=OFF \
-DTpetra_INST_COMPLEX_DOUBLE=OFF \
-DTpetra_INST_INT_LONG_LONG:BOOL=OFF \
-DTPL_ENABLE_Boost:BOOL=OFF \
-DTPL_ENABLE_SuperLU:BOOL=OFF \
-DTrilinos_ALLOW_NO_PACKAGES:BOOL=OFF \
-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF  \
-DTrilinos_ENABLE_AztecOO:BOOL=OFF \
-DTrilinos_ENABLE_Epetra:BOOL=OFF \
-DTrilinos_ENABLE_EpetraExt:BOOL=OFF \
-DTrilinos_ENABLE_Fortran:BOOL=OFF \
-DTrilinos_ENABLE_Gtest:BOOL=OFF \
-DTrilinos_ENABLE_MueLu:BOOL=OFF \
-DTrilinos_ENABLE_NOX:BOOL=OFF  \
-DTrilinos_ENABLE_STKUnit_test_utils:BOOL=OFF \
-DTrilinos_ENABLE_STKUnit_tests:BOOL=OFF \
-DTrilinos_ENABLE_TESTS:BOOL=OFF  \
-DKokkos_BINARY_DIR=/user-environment/env/default/bin -DKokkos_CUDA_DIR=/user-environment/env/default/

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
	-DUTOPIA_INSTALL_YAML_CPP=ON \
	-DPETSC_DIR=$INSTALL_DIR/petsc \
	-DTrilinos_DIR=$INSTALL_DIR/Trilinos/lib64/cmake/Trilinos

make yaml-cpp

cmake .. 
make -j72 && make install
