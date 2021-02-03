#!/bin/bash
# utopia_install_spack.sh

DEFAULT_SPACK_DIR=$PWD/spack
SPACK_DIR=${1:-$DEFAULT_SPACK_DIR}
git clone https://github.com/spack/spack $SPACK_DIR
. $SPACK_DIR/share/spack/setup-env.sh && export PATH=$SPACK_DIR/bin:$PATH

spack install trilinos +tpetra +openmp +muelu +amesos2 +ifpack2 +ml #+adios2
spack install doxygen python swig

# Install pytorch
spack install anaconda3
spack load anaconda3
conda install --yes  pytorch torchvision torchaudio cpuonly -c pytorch

# spack install petsc +superlu-dist +mumps +metis +exodusii