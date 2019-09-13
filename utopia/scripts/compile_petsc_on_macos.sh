#!/bin/bash

./configure --prefix=$INSTALL_DIR/petsc --with-mpi=1 --download-superlu_dist=yes --download-superlu=yes --download-mumps=yes --download-scalapack=yes --download-hypre=yes --with-mpi-dir=/opt/local --with-cxx-dialect=C++11 --with-debugging=0 
#--download-slepc=yes