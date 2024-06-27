#!/bin/bash                                                                                                                                                                                           

# source /opt/cray/pe/cpe/23.12/restore_lmod_system_defaults.sh                                                                                                                                       

module load cray
module load PrgEnv-gnu
module load CMake
module load cray-mpich
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel

export CXX=CC
export CC=cc

# export FC=ftn
# export F90=ftn
# export F77=ftn
