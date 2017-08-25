export MOOSE_DIR=/apps/daint/UES/jenkins/6.0.UP02/gpu/easybuild/software/moose/d2193ae-CrayGNU-2016.11
export LIBMESH_DIR=$MOOSE_DIR/libmesh/installed

#utopia modules
module unload PrgEnv-cray
module load PrgEnv-gnu 
module load nano
module load CMake
module load slurm
module load daint-mc
module load cray-petsc

#utopia_fe modules
module load cray-netcdf 

