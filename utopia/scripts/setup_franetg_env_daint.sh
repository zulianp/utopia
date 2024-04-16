#!/bin/bash

# module load daint-gpu

# TODO: Add if statement.
# module switch PrgEnv-cray PrgEnv-gnu
# module load CMake/3.26.5
# module switch gcc/11.2.0 gcc/8.3.0
# module load cray-mpich
# module load cray-hdf5-parallel
# module load cray-netcdf-hdf5parallel
# module load cudatoolkit
# module load cuda_memtest


source $APPS/UES/anfink/gpu/environment


# Temporary.
export STAGE_DIR=$SCRATCH/code
export INSTALL_DIR=$STAGE_DIR/installations
export UTOPIA_DIR=$INSTALL_DIR/utopia_franetg
export MOONOLITH_DIR=$INSTALL_DIR/par_moonolith
export UTOPIA_FE_EXEC=$INSTALL_DIR/utopia_fe_franetg/bin/utopia_fe_exec
export TRILINOS_DIR=$STAGE_DIR/utopia/external/Trilinos
export ADIOS_DIR=$INSTALL_DIR/adios2
export SuperLU_DIR=$INSTALL_DIR/superlu
export CRAYPE_LINK_TYPE=dynamic
export PNETCDF_DIR=/opt/cray/pe/parallel-netcdf/1.12.1.4/gnu/8.2
export HDF5_DIR=$HDF5_ROOT