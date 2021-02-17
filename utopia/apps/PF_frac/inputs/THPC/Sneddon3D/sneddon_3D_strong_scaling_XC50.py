import os
import subprocess

# uncomment the lines for the small/big experiment that you want to submit jobs into the queue
# inputFile is the full path to the file containing the experiments parameters
# basedir is the starting directory where each experiment setup will be saved (can be relative or absolute)

# num_nodes = (4, 5, 6, 7, 8, 12, 16, 20, 24, 28, 32)
num_nodes = (1, 2, 4, 8)
inputFile = '$SCRATCH/Sneddon_small.json'
basedir = 'Sneddon3D_XC50_small'
# num_nodes = (40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256)
# inputFile = '$SCRATCH/Sneddon3D.json'
# basedir = 'Sneddon3D_XC50_large'


for nodes in num_nodes:
    dirname = f'{basedir}/{nodes:03}nodes'
    jobfile = 'job.sh'
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(os.path.join(dirname, jobfile), 'w') as scriptfile:
        scriptfile.write(f"""#!/bin/bash
#SBATCH --job-name=SS_Sneddon3D_XC50_small
#SBATCH --output=output_%j.out
#SBATCH --time=00:50:00
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --constraint=gpu
#SBATCH --account=sm43

export OMP_PROC_BIND=false
export OMP_NUM_THREADS=1
source $APPS/UES/anfink/cpu/environment
srun $SCRATCH/utopia_exec @file {inputFile}""")

    subprocess.call(['sbatch', jobfile], cwd=dirname)
