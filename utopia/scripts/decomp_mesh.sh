n_procs=$1

$TRILINOS_DIR/bin/decomp -p $n_procs mesh_matrix.e
$TRILINOS_DIR/bin/decomp -p $n_procs mesh_fracture.e