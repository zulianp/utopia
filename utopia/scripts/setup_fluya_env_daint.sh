#!/bin/bash

module load daint-mc
module unload PrgEnv-cray
module load PrgEnv-gnu
module load cray-mpich
module load cray-hdf5
module load cray-netcdf 
module load cray-parallel-netcdf
module load CMake/3.26.5

export CXX=CC
export CC=cc
export FC=ftn
export F90=ftn
export F77=ftn
