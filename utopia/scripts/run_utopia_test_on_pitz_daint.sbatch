#!/bin/bash -l
#SBATCH --job-name=UtopiaTest
#SBATCH --mail-user=patrick.zulian@usi.ch
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --partition=debug
#SBATCH --constraint=mc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun ../bin/utopia_exec -on_error_abort -verbose -bench
