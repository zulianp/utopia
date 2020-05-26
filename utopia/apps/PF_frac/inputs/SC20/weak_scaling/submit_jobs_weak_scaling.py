import json
import math
import os
import subprocess

# inputFile is the full path to the file containing the experiments parameters
# basedir is the starting directory where each experiment setup will be saved (can be relative or absolute)

num_nodes = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
inputFile = os.path.join(os.environ['SCRATCH'], 'build/paper_sc/anfink/inputs_sc20/frac_net3D_weak.json')
basedir = f'weak_scaling'

inputFileData = json.load(open(inputFile))

for nodes in num_nodes:
    dirname = f'{basedir}/{nodes:04}nodes'
    jobfile = 'job.sh'
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(os.path.join(dirname, jobfile), 'w') as scriptfile:
        scriptfile.write(f"""#!/bin/bash
#SBATCH --job-name=utopia_weak_scaling
#SBATCH --output=output_%j.out
#SBATCH --time=00:30:00
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --constraint=gpu

export OMP_PROC_BIND=false
export OMP_NUM_THREADS=1
source $APPS/UES/anfink/cpu_20200218_Release/environment
srun $SCRATCH/build/cpu/utopia_edsl_refactor/build/utopia_exec @file frac_net3D_weak.json""")

    domainSize = math.ceil(math.pow(1000*nodes, 1/3))
    inputFileData['nx'] = domainSize
    inputFileData['ny'] = domainSize
    inputFileData['nz'] = domainSize
    json.dump(inputFileData, open(os.path.join(dirname, 'frac_net3D_weak.json'), 'w'))

    subprocess.call(['sbatch', jobfile], cwd=dirname)

