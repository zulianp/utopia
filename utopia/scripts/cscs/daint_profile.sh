#!/bin/bash
  
if [[ "$-" == *i* ]]; then
       bind '"\e[A":history-search-backward'
       bind '"\e[B":history-search-forward'
fi

module switch PrgEnv-cray PrgEnv-gnu
module load CMake
module load git
module load cray-hdf5-parallel
module load cudatoolkit
module load cuda_memtest
module load daint-gpu

export CRAYPE_LINK_TYPE=dynamic
export CXX=$TRILINOS_DIR/bin/nvcc_wrapper
export SRCDIR=$HOME/utopia/
export UTOPIA_DIR=$HOME/utopia/build/
source /scratch/snx1600/anfink/shared_pasc/20180919_eurohack18/environment

export CXX=$TRILINOS_DIR/bin/nvcc_wrapper
export CUDA_LAUNCH_BLOCKING=1
