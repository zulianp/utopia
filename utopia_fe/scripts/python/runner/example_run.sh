#!/bin/bash

python3 utopia_fe.py \
    --solid=./solid.e  \
    --fluid=./fluid.e  \
    --workspace=./workspace  \
    --young_modulus=9.174e6  \
    --poisson_ratio=0.4  \
    --barrier_parameter=3e-9  \
    --min_barrier_parameter=3e-9  \
    --barrier_parameter_shrinking_factor=0.3 \
    --mass=998 \
    --containement_iterations=40 \
    -a \
    -p 4


