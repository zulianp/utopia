#!/bin/bash
# utopia_install_spack.sh

DEFAULT_SPACK_DIR=$PWD/spack
SPACK_DIR=${1:-$DEFAULT_SPACK_DIR}

git clone https://github.com/spack/spack $SPACK_DIR && \
    . $SPACK_DIR/share/spack/setup-env.sh && export PATH=$SPACK_DIR/bin:$PATH && \
    spack install cmake doxygen clang-format gfortran python swig &&  \
    spack install anaconda3 &&  \
    spack install openblas  &&  \
    spack install lapack    &&  \
    spack load anaconda3    &&  \
    conda install --yes  pytorch torchvision torchaudio cpuonly -c pytorch

# spack install petsc +superlu-dist +mumps +metis +exodusii
# spack install trilinos +tpetra +openmp +muelu +amesos2 +adios2 +ifpack2 +ml