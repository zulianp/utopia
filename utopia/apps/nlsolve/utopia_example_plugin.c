#include "utopia_plugin_Function.h"

#include <stdio.h>

// How to create a plugin (inside the build folder)
// mpicc -c ../apps/nlsolve/utopia_example_plugin.c -I ../backend/plugin
// mpicc -shared utopia_example_plugin.o -o utopia_example_plugin.dylib

int UTOPIA_PLUGIN_EXPORT utopia_load_plugin(const char *path) {
    printf("utopia_load_plugin\n");
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_init(MPI_Comm comm, plugin_Function_t *info) {
    printf("utopia_plugin_Function_init\n");
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_crs_graph(ptrdiff_t *nlocal,
                                                                 ptrdiff_t *nglobal,
                                                                 ptrdiff_t *nnz,
                                                                 plugin_idx_t **rowptr,
                                                                 plugin_idx_t **colidx) {
    printf("utopia_plugin_Function_create_crs_graph\n");
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_vector(ptrdiff_t *nlocal,
                                                              ptrdiff_t *nglobal,
                                                              plugin_scalar_t **values) {
    printf("utopia_plugin_Function_create_vector\n");
    return UTOPIA_PLUGIN_FAILURE;
}

// Optimization function
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_value(const plugin_scalar_t *x, plugin_scalar_t *const out) {
    printf("utopia_plugin_Function_value\n");
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_gradient(const plugin_scalar_t *const x, plugin_scalar_t *const out) {
    printf("utopia_plugin_Function_gradient\n");
    return UTOPIA_PLUGIN_FAILURE;
}

// We might want to have also other formats here
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_hessian_crs(const plugin_scalar_t *const x,
                                                            const plugin_idx_t *const rowptr,
                                                            const plugin_idx_t *const colidx,
                                                            plugin_scalar_t *const values) {
    printf("utopia_plugin_Function_hessian_crs\n");
    return UTOPIA_PLUGIN_FAILURE;
}

// Operator
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply(const plugin_scalar_t *const x, plugin_scalar_t *const out) {
    printf("utopia_plugin_Function_apply\n");
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply_constraints(const plugin_scalar_t *const x) {
    printf("utopia_plugin_Function_apply_constraints\n");
    return UTOPIA_PLUGIN_FAILURE;
}
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply_zero_constraints(const plugin_scalar_t *const x) {
    printf("utopia_plugin_Function_apply_zero_constraints\n");
    return UTOPIA_PLUGIN_FAILURE;
}
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_copy_constrained_dofs(const plugin_scalar_t *const src,
                                                                      plugin_scalar_t *const dest) {
    printf("utopia_plugin_Function_copy_constrained_dofs\n");
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy(plugin_Function_t *info) {
    printf("utopia_plugin_Function_destroy\n");
    return UTOPIA_PLUGIN_FAILURE;
}
