#!/bin/bash

# module switch PrgEnv-cray PrgEnv-gnu
# module load CMake
# module load cray-mpich
# module load cray-hdf5-parallel
# module load cray-netcdf-hdf5parallel
# export CC=cc
# export CXX=CC
# export FC=ftn
# export F90=ftn
# export F77=ftn

source $APPS/UES/anfink/cpu/environment

# We compile the new versions!
unset TRILINOS_DIR
unset PETSC_DIR
unset LIBMESH_DIR
unset MOONOLITH_DIR
unset UTOPIA_DIR