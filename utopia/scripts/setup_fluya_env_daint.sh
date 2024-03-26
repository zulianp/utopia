#!/bin/bash

module load daint-mc
module switch PrgEnv-cray PrgEnv-gnu
module load CMake
module load gcc
module load cray-mpich
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel




# Temporary.
export STAGE_DIR=$SCRATCH/code
export INSTALL_DIR=$STAGE_DIR/installations
export UTOPIA_DIR=$INSTALL_DIR/utopia_fluya
export MOONOLITH_DIR=$INSTALL_DIR/par_moonolith
export UTOPIA_FE_EXEC=$INSTALL_DIR/utopia_fluya/utopia_fe_exec
export TRILINOS_DIR=$STAGE_DIR/utopia/external/Trilinos