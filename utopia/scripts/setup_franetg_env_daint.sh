#!/bin/bash

# Setup environment for franetg.
source /apps/daint/UES/anfink/gpu_2023-04-17_Release/environment

# Temporary.
export STAGE_DIR=$SCRATCH/code
export INSTALL_DIR=$STAGE_DIR/installations
export METIS_DIR=$INSTALL_DIR/metis
export UTOPIA_DIR=$INSTALL_DIR/utopia_franetg_gpu
export UTOPIA_FE_EXEC=$INSTALL_DIR/utopia_fe_franetg_gpu/bin/utopia_fe_exec
export ADIOS_DIR=$INSTALL_DIR/adios2_gpu
export SuperLU_DIR=$INSTALL_DIR/superlu
export CRAYPE_LINK_TYPE=dynamic
export HDF5_DIR=$HDF5_ROOT