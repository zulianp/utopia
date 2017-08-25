#in order to lunch this connect to the cluster with ssh -X username@ela.cscs.ch and ssh -X daint

ALLINEA_MPI_INIT_PENDING=1 
ALLINEA_MPI_INIT=MPI_Init_thread

current=$PWD
cp ../bin/utopia_exec $SCRATCH/
cp run_utopia_test_on_pitz_daint.sbatch $SCRATCH/
cd $SCRATCH

sbatch run_utopia_test_on_pitz_daint.sbatch


echo "use cd "$current" to go to script folder"