#!/bin/bash


# ./configure --prefix=$INSTALL_DIR/petsc --with-packages-search-path=/opt/local  --with-mpi=1 --with-hypre=1 --with-superlu_dist=1 --with-superlu=1   --download-mumps=yes  --download-scalapack=yes --download-libmesh=yes 


./configure --prefix=$INSTALL_DIR/petsc --with-mpi=1 --with-hypre=1 --with-superlu_dist=1 --with-superlu=1   --download-mumps=yes  --download-scalapack=yes --with-hypre-dir=/opt/local --with-superlu_dist-dir=/opt/local --with-superlu-dir=/opt/local   --download-mumps=yes
