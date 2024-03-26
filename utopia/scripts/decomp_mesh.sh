n_procs=$1

srun $TRILINOS_DIR/bin/decomp -p $n_procs mesh_matrix.e
srun $TRILINOS_DIR/bin/decomp -p $n_procs mesh_fracture.e