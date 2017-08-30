module load daint-mc
# module unload daint-mc
# module load daint-gpu

script_dir=$PWD
cd $SCRATCH && mkdir bin && cd bin

cp $script_dir/../bin/utopia_fe_exec $SCRATCH/
cp $script_dir/../data $SCRATCH/

sbatch $script_dir/utopia_fe_experiment.sbatch

squeue -u zulianp

echo "-----------------------------------------"
echo "use cd $script_dir to go to script folder"
echo "-----------------------------------------"
