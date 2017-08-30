#!/bin/bash
ntasks=40
ntasks_per_node=10
partition=slim
duration=1:00:00

sbatch  --ntasks=${ntasks}  --ntasks-per-node=${ntasks_per_node} --partition=${partition} --time=${duration} utopia_fe_experiment_on_ics_hpc.job
squeue -u zulian
