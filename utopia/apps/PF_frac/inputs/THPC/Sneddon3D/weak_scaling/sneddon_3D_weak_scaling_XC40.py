import os
import subprocess

num_nodes = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
tasks_per_node=36

inputFile = '$SCRATCH/sneddon_3D_weak_scaling.json'
basedir = 'Sneddon3D_weak'

for nodes in num_nodes:
    dirname = f'{basedir}/{nodes:03}nodes'
    jobfile = 'job.sh'
    domainSize = math.ceil(math.pow(1000*nodes, 1/3))

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(os.path.join(dirname, jobfile), 'w') as scriptfile:
        scriptfile.write(f"""#!/bin/bash
#SBATCH --job-name=WS_Sneddon3D_weak
#SBATCH --output=output_%j.out
#SBATCH --time=00:10:00
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --account=u2

export OMP_PROC_BIND=false
export OMP_NUM_THREADS=1
source $APPS/UES/anfink/cpu/environment
srun $SCRATCH/utopia_exec @file {inputFile} -nx {domainSize} -ny {domainSize} -nz {domainSize}""")

    subprocess.call(['sbatch', jobfile], cwd=dirname)
