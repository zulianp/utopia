import os
import subprocess

# uncomment the lines for the small/big experiment that you want to submit jobs into the queue
# inputFile is the full path to the file containing the experiments parameters
# basedir is the starting directory where each experiment setup will be saved (can be relative or absolute)

#num_nodes = (4, 5, 6, 7, 8, 12, 16, 20, 24, 28, 32)
#inputFile = '$SCRATCH/build/paper_sc/anfink/inputs_sc20/frac_net3D_25x25x25.json'
#basedir = '25x25x25_5'
num_nodes = (40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256)
inputFile = '$SCRATCH/build/paper_sc/anfink/inputs_sc20/frac_net3D_50x50x50.json'
basedir = '50x50x50'

for nodes in num_nodes:
    dirname = f'{basedir}/{nodes:03}nodes'
    jobfile = 'job.sh'
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(os.path.join(dirname, jobfile), 'w') as scriptfile:
        scriptfile.write(f"""#!/bin/bash
#SBATCH --job-name=utopia_strong_scaling
#SBATCH --output=output_%j.out
#SBATCH --time=02:30:00
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --constraint=gpu

export OMP_PROC_BIND=false
export OMP_NUM_THREADS=1
source $APPS/UES/anfink/cpu_20200218_Release/environment
srun $SCRATCH/build/cpu/utopia_edsl_refactor/build/utopia_exec @file {inputFile}""")

    subprocess.call(['sbatch', jobfile], cwd=dirname)

