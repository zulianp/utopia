#!/bin/bash
# utopia_env.sh

DEFAULT_SPACK_DIR=$PWD/spack
SPACK_DIR=${1:-$DEFAULT_SPACK_DIR}
. $SPACK_DIR/share/spack/setup-env.sh

spack load trilinos
spack load anaconda3
spack load cmake doxygen clang-format gfortran python swig 