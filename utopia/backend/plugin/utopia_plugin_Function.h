#ifndef UTOPIA_PLUGIN_FUNCTION_H
#define UTOPIA_PLUGIN_FUNCTION_H

#include <mpi.h>

typedef struct {
    MPI_Comm comm;
    void *user_data;
} plugin_Function_t;

typedef double plugin_scalar_t;
typedef int plugin_idx_t;

#define UTOPIA_PLUGIN_EXPORT __attribute__((visibility("default")))

enum UtopiaPluginCode { UTOPIA_PLUGIN_SUCCESS = 0, UTOPIA_PLUGIN_NAN, UTOPIA_PLUGIN_FAILURE };

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy_array(const plugin_Function_t *info, void *ptr);
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_array(const plugin_Function_t *info, size_t size, void **ptr);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_init(plugin_Function_t *info);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_crs_graph(const plugin_Function_t *info,
                                                                 ptrdiff_t *nlocal,
                                                                 ptrdiff_t *nglobal,
                                                                 ptrdiff_t *nnz,
                                                                 plugin_idx_t **rowptr,
                                                                 plugin_idx_t **colidx);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_create_vector(const plugin_Function_t *info,
                                                              ptrdiff_t *nlocal,
                                                              ptrdiff_t *nglobal,
                                                              plugin_scalar_t **values);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy_vector(const plugin_Function_t *info, plugin_scalar_t *values);

// Optimization function
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_value(const plugin_Function_t *info,
                                                      const plugin_scalar_t *x,
                                                      plugin_scalar_t *const out);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_gradient(const plugin_Function_t *info,
                                                         const plugin_scalar_t *const x,
                                                         plugin_scalar_t *const out);

// We might want to have also other formats here
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_hessian_crs(const plugin_Function_t *info,
                                                            const plugin_scalar_t *const x,
                                                            const plugin_idx_t *const rowptr,
                                                            const plugin_idx_t *const colidx,
                                                            plugin_scalar_t *const values);

// Operator
int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply(const plugin_Function_t *info,
                                                      const plugin_scalar_t *const x,
                                                      plugin_scalar_t *const out);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply_constraints(const plugin_Function_t *info,
                                                                  plugin_scalar_t *const x);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_apply_zero_constraints(const plugin_Function_t *info,
                                                                       plugin_scalar_t *const x);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_copy_constrained_dofs(const plugin_Function_t *info,
                                                                      const plugin_scalar_t *const src,
                                                                      plugin_scalar_t *const dest);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_report_solution(const plugin_Function_t *info,
                                                                const plugin_scalar_t *const x);

int UTOPIA_PLUGIN_EXPORT utopia_plugin_Function_destroy(plugin_Function_t *info);

#endif  // UTOPIA_PLUGIN_FUNCTION_H
