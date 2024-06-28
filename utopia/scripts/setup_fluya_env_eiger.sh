#!/bin/bash                                                                                                                                                                                           

module load cray
module unload PrgEnv-cray
module load PrgEnv-gnu
module load cray-mpich
module load cray-hdf5
module load cray-netcdf 
module load cray-parallel-netcdf
module load cray-libsci
module load CMake

export CXX=CC
export CC=cc
export FC=ftn
export F90=ftn
export F77=ftn


# module load cray-hdf5-parallel
# module load cray-netcdf-hdf5parallel