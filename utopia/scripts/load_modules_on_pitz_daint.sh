export MOOSE_DIR=/apps/daint/UES/jenkins/6.0.UP02/gpu/easybuild/software/moose/d2193ae-CrayGNU-2016.11
export LIBMESH_DIR=$MOOSE_DIR/libmesh/installed

#a decent text editor
module load nano

#utopia modules
module unload PrgEnv-cray
module load PrgEnv-gnu 
module load daint-mc
module load CMake
module load slurm
module load cray-petsc

module load cray-libsci
module load cray-libsci_acc
module load craype-accel-nvidia35
module load cray-tpsl
module load Boost
module load GSL
module load Visit
module load ddt

#utopia_fe modules
# module load cray-netcdf #this breaks the linking with libmesh

# export compiler wrapper paths
# export CC=cc
# export CXX=CC
# export FC=ftn
# export F90=ftn
# export F77=ftn

# make dynamic linking the default
export CRAYPE_LINK_TYPE=dynamic
export CRAY_ADD_RPATH=yes
export CRAY_CPU_TARGET=x86-64
