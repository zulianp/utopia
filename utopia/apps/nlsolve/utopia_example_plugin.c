#include "utopia_plugin_Function.h"

#include <stdio.h>
#include <stdlib.h>

// How to create a plugin (inside the build folder)
// mpicc -c ../apps/nlsolve/utopia_example_plugin.c -I ../backend/plugin
// mpicc -shared utopia_example_plugin.o -o utopia_example_plugin.dylib

int UTOPIA_PLUGIN_EXPORT utopia_load_plugin(const char *path) {
    printf("utopia_load_plugin\n");
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_init(MPI_Comm comm, plugin_Function_t *info) {
    printf("utopia_plugin_Function_init\n");
    int size;
    MPI_Comm_size(comm, &size);
    if (size != 1)
        return UTOPIA_PLUGIN_FAILURE;
    else
        return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_crs_graph(ptrdiff_t *nlocal,
                                                                 ptrdiff_t *nglobal,
                                                                 ptrdiff_t *nnz,
                                                                 plugin_idx_t **rowptr,
                                                                 plugin_idx_t **colidx) {
    printf("utopia_plugin_Function_create_crs_graph\n");
    *nlocal = 2;
    *nglobal = 2;
    *nnz = 2;
    // TODO
    return UTOPIA_PLUGIN_FAILURE;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_vector(ptrdiff_t *nlocal,
                                                              ptrdiff_t *nglobal,
                                                              plugin_scalar_t **values) {
    *values = (plugin_scalar_t *)malloc(2 * sizeof(plugin_scalar_t));
    (*values)[0] = -1;
    (*values)[1] = -1;
    *nlocal = 2;
    *nglobal = 2;
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy_vector(plugin_scalar_t *values) {
    free(values);
    return UTOPIA_PLUGIN_SUCCESS;
}

// Optimization function
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_value(const plugin_scalar_t *x, plugin_scalar_t *const out) {
    *out = (x[0] * x[0] + x[1] * x[1]);
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_gradient(const plugin_scalar_t *const x, plugin_scalar_t *const out) {
    out[0] = 2 * x[0];
    out[1] = 2 * x[1];
    return UTOPIA_PLUGIN_SUCCESS;
}

// We might want to have also other formats here
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_hessian_crs(const plugin_scalar_t *const x,
                                                            const plugin_idx_t *const rowptr,
                                                            const plugin_idx_t *const colidx,
                                                            plugin_scalar_t *const values) {
    printf("utopia_plugin_Function_hessian_crs\n");

    values[0] = 2;
    values[1] = 2;

    return UTOPIA_PLUGIN_SUCCESS;
}

// Operator
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply(const plugin_scalar_t *const x, plugin_scalar_t *const out) {
    printf("utopia_plugin_Function_apply\n");
    out[0] = 2 * x[0];
    out[1] = 2 * x[1];
    return UTOPIA_PLUGIN_SUCCESS;
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
