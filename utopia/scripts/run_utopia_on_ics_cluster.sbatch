#!/bin/bash -l
#SBATCH --job-name=UtopiaTest
#SBATCH --mail-type=ALL
#SBATCH --mail-user=patrick.zulian@usi.ch
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1

srun /home/zulian/utopia/utopia/bin/utopia_exec -on_error_abort --verbose
