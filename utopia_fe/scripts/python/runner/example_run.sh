#!/bin/bash

python utopia_fe.py \
    --solid=solid_mesh.e  \
    --fluid=fluid_mesh.e  \
    --workspace=./my_output_dir  \
    --young_modulus=9.174e6  \
    --poisson_ratio=0.4  \
    --barrier_parameter=1e-12  \
    --min_barrier_parameter=1e-8  \
    --barrier_parameter_shrinking_factor=0.3 \
    --mass=998 \
    --containement_iterations=40 \
    -g # With -g it generates YAML files.
       # Use -a for running all stages or -r <stage_num>
       # for running a particular stage.

