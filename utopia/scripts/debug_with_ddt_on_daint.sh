# For lunching this script connect to the cluster with ssh -X username@ela.cscs.ch and ssh -X daint
# Eventually disconnecting and reconnecting with the cluster might fix problems arising with X server

module load daint-mc
module load ddt
export ALLINEA_MPI_INIT_PENDING=1 
export ALLINEA_MPI_INIT=MPI_Init_thread

script_dir=$PWD
cp ../bin/utopia_exec $SCRATCH/
cd $SCRATCH

salloc -n 8 --time=00:30:00 -C mc --partition=debug
ddt --processes=8 ./utopia_exec 

echo "use cd $script_dir to go to script folder"
