#!/bin/bash
  
if [[ "$-" == *i* ]]; then
       bind '"\e[A":history-search-backward'
       bind '"\e[B":history-search-forward'
fi

module switch PrgEnv-cray PrgEnv-gnu
module load CMake/3.8.1
module load git/2.11.0
module load cray-hdf5-parallel
module load cudatoolkit
module load cuda_memtest/1.2.3
module load daint-gpu
export CXX=$TRILINOS_DIR/bin/nvcc_wrapper
export SRCDIR=/users/hck40/utopia/
export UTOPIA_DIR=/users/hck40/utopia/build/
source /scratch/snx1600/anfink/shared_pasc/20180919_eurohack18/environment

export CXX=$TRILINOS_DIR/bin/nvcc_wrapper
export CUDA_LAUNCH_BLOCKING=1
