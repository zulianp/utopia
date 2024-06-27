#!/bin/bash

module load daint-mc
module unload PrgEnv-cray
module load PrgEnv-gnu
module load cray-hdf5 cray-netcdf
module load cray-mpich
module load CMake/3.26.5

export CXX=CC
export CC=cc
export FC=ftn
export F90=ftn
export F77=ftn
