#!/usr/bin/env bash

# This example requires the function plugin from sfem

set -e
set -x

export SFEM_MESH_DIR=/Users/patrickzulian/Desktop/code/utopia/utopia_fe/data/hydros/mesh-multi-outlet-better
export SFEM_DIRICHLET_NODES=$SFEM_MESH_DIR/zd.raw 
export SFEM_NEUMAN_FACES=$SFEM_MESH_DIR/on.raw  
PLUGIN=/Users/patrickzulian/Desktop/code/sfem/isolver_sfem.dylib 

time ./utopia_exec -app nlsolve -path $PLUGIN --verbose -solver_type Newton -max_it 10 -rtol 1e-10 -atol 1e-16 -stol 1e-16
