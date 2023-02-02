#include "utopia_plugin_Function.h"

#include <stdio.h>
#include <stdlib.h>

// How to create a plugin (inside the build folder)
// mpicc -c ../apps/nlsolve/utopia_example_plugin.c -I ../backend/plugin
// mpicc -shared utopia_example_plugin.o -o utopia_example_plugin.dylib

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy_array(const plugin_Function_t *info, void *ptr) {
    free(ptr);
    return UTOPIA_PLUGIN_SUCCESS;
}
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_array(const plugin_Function_t *info, size_t size, void **ptr) {
    *ptr = malloc(size);
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_init(plugin_Function_t *info) {
    int size;
    MPI_Comm_size(info->comm, &size);
    if (size != 1) {
        // Only serial runs for this plugin!
        return UTOPIA_PLUGIN_FAILURE;
    } else {
        // Eventuall initialize info->user_data here
        return UTOPIA_PLUGIN_SUCCESS;
    }
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_crs_graph(const plugin_Function_t *info,
                                                                 ptrdiff_t *nlocal,
                                                                 ptrdiff_t *nglobal,
                                                                 ptrdiff_t *nnz,
                                                                 plugin_idx_t **rowptr,
                                                                 plugin_idx_t **colidx) {
    *nlocal = 2;
    *nglobal = 2;
    *nnz = 2;
    // TODO
    *rowptr = (plugin_idx_t *)malloc((*nlocal + 1) * sizeof(plugin_idx_t));
    (*rowptr)[0] = 0;
    (*rowptr)[1] = 1;
    (*rowptr)[2] = 2;

    *colidx = (plugin_idx_t *)malloc((*nnz) * sizeof(plugin_idx_t));

    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_vector(const plugin_Function_t *info,
                                                              ptrdiff_t *nlocal,
                                                              ptrdiff_t *nglobal,
                                                              plugin_scalar_t **values) {
    *values = (plugin_scalar_t *)malloc(2 * sizeof(plugin_scalar_t));
    (*values)[0] = -1;
    (*values)[1] = -1;
    *nlocal = 2;
    *nglobal = 2;
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy_vector(const plugin_Function_t *info, plugin_scalar_t *values) {
    free(values);
    return UTOPIA_PLUGIN_SUCCESS;
}

// Optimization function
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_value(const plugin_Function_t *info,
                                                      const plugin_scalar_t *x,
                                                      plugin_scalar_t *const out) {
    *out = (x[0] * x[0] + x[1] * x[1]);
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_gradient(const plugin_Function_t *info,
                                                         const plugin_scalar_t *const x,
                                                         plugin_scalar_t *const out) {
    out[0] = 2 * x[0];
    out[1] = 2 * x[1];
    return UTOPIA_PLUGIN_SUCCESS;
}

// We might want to have also other formats here
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_hessian_crs(const plugin_Function_t *info,
                                                            const plugin_scalar_t *const x,
                                                            const plugin_idx_t *const rowptr,
                                                            const plugin_idx_t *const colidx,
                                                            plugin_scalar_t *const values) {
    values[0] = 2;
    values[1] = 2;
    return UTOPIA_PLUGIN_SUCCESS;
}

// Operator
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply(const plugin_Function_t *info,
                                                      const plugin_scalar_t *const x,
                                                      const plugin_scalar_t *const h,
                                                      plugin_scalar_t *const out) {
    out[0] = 2 * h[0];
    out[1] = 2 * h[1];
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply_constraints(const plugin_Function_t *info,
                                                                  plugin_scalar_t *const x) {
    // No constraints for this example
    return UTOPIA_PLUGIN_SUCCESS;
}
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply_zero_constraints(const plugin_Function_t *info,
                                                                       plugin_scalar_t *const x) {
    // No constraints for this example
    return UTOPIA_PLUGIN_SUCCESS;
}
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_copy_constrained_dofs(const plugin_Function_t *info,
                                                                      const plugin_scalar_t *const src,
                                                                      plugin_scalar_t *const dest) {
    // No constraints for this example
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_report_solution(const plugin_Function_t *info,
                                                                const plugin_scalar_t *const x) {
    printf("%s: Result %g %g\n", __FILE__, x[0], x[1]);
    return UTOPIA_PLUGIN_SUCCESS;
}

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy(plugin_Function_t *info) {
    // No user-data for this example
    return UTOPIA_PLUGIN_SUCCESS;
}
