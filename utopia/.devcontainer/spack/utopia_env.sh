#!/bin/bash
# utopia_env.sh


DEFAULT_SPACK_DIR=$PWD/spack
SPACK_DIR=${1:-$DEFAULT_SPACK_DIR}

. $SPACK_DIR/share/spack/setup-env.sh

spack load cmake doxygen clang-format gfortran python swig 
spack load anaconda3
spack load openblas
spack load lapack
# spack load petsc
# spack load trilinos