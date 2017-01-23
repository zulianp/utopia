#/bin/bash

#tested with petsc compiled with ./configure --with-cxx=0 --with-fc=0 --with-mpi-dir=$MPICH_DIR  --download-f2cblaslapack=1 --with-shared-libraries=0 --with-pic=1

export PETSC_ROOT=$MOONOLITH_ROOT/external/petsc
export PETSC_DIR=$PETSC_ROOT
export PETSC_INCLUDE_CONF=$PETSC_ROOT/arch-linux2-c-debug/include
export PETSC_ARCH=arch-linux2-c-debug

cd ../bin && cmake .. -DPETSC_INCLUDES=$PETSC_DIR/include -DPETSC_INCLUDE_CONF=$PETSC_INCLUDE_CONF -DPETSC_LIBRARIES=$PETSC_DIR/lib/libpetsc.dylib -DPETSC_LIBRARY_SINGLE=$PETSC_DIR/lib/libpetsc.dylib -DPETSC_EXECUTABLE_RUNS=ON -DPETSC_COMPILER=/opt/cray/craype/2.4.0/bin/CC -DWITH_PETSC_CRAY=ON && make -j4
cd ../scripts