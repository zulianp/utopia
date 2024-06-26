#!/bin/bash -l
set -e
set -x

#SBATCH --job-name=large_tutorials_utopia_fe
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --account=u0
#SBATCH --constraint=gpu

# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# srun ../bin/utopia_exec -on_error_abort -verbose -bench

if [[ -z "$UTOPIA_FE_EXEC" ]]
then
        echo "define UTOPIA_FE_EXEC"
        exit 1
fi

if [[ -z "$TRILINOS_DIR" ]]
then
        echo "define TRILINOS_DIR"
        exit 1
fi


./decomp_mesh.sh 12

srun $UTOPIA_FE_EXEC @file ../flow.yaml
# srun $UTOPIA_FE_EXEC @file ../transport.yaml
