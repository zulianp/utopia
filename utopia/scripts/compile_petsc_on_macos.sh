#!/bin/bash

./configure --prefix=$INSTALL_DIR/petsc --with-mpi=1 --with-hypre=1 --download-superlu_dist=1 --with-superlu=1 --download-mumps=yes --download-scalapack=yes --with-hypre-dir=/opt/local --with-superlu-dir=/opt/local --download-mumps=yes --with-mpi-dir=/opt/local --with-cxx-dialect=C++11 --with-debugging=0 
#--download-slepc=yes