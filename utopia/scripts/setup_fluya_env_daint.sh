#!/bin/bash

# Setup environment for daint fluya.
source /apps/daint/UES/anfink/gpu_2023-04-17_Release/environment

# Temporary.
export STAGE_DIR=$SCRATCH/code
export INSTALL_DIR=$STAGE_DIR/installations
export METIS_DIR=$INSTALL_DIR/metis
export UTOPIA_DIR=$INSTALL_DIR/utopia_fluya_gpu
export UTOPIA_FE_EXEC=$INSTALL_DIR/utopia_fe_fluya_gpu/bin/utopia_fe_exec
export SuperLU_DIR=$INSTALL_DIR/superlu
export CRAYPE_LINK_TYPE=dynamic
export HDF5_DIR=$HDF5_ROOT


# module load daint-mc

# # TODO: Add if statement.
# module switch PrgEnv-cray PrgEnv-gnu
# module load CMake
# module load gcc
# module load cray-mpich
# module load cray-hdf5-parallel
# module load cray-netcdf-hdf5parallel


# # Temporary.
# export STAGE_DIR=$SCRATCH/code
# export INSTALL_DIR=$STAGE_DIR/installations
# export UTOPIA_DIR=$INSTALL_DIR/utopia_fluya
# export MOONOLITH_DIR=$INSTALL_DIR/par_moonolith
# export UTOPIA_FE_EXEC=$INSTALL_DIR/utopia_fe_fluya/bin/utopia_fe_exec
# export TRILINOS_DIR=$STAGE_DIR/utopia/external/Trilinos
# export KOKKO_DIR=$TRILINOS_DIR/lib/cmake/Kokkos