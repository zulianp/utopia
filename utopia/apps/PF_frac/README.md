# Phase-field simulations

Phase field simulations now support checkpoints and simulation restart.
Add these three lines to your `.yaml` file.
```yaml
restart: true					# Enable the checkpoint/restart functionality
checkpoint_each: 10				# Create a checkpoint each 10 steps (for slow simulations it is ok also to use 1)
checkpoint_path: checkpoint		# Folder where checkpoints are saved. Delete folder if you do not want to restart from these checkpoints
```

**The functionaly has been tested only if the number of MPI processes accross runs is kept the same.**
