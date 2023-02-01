#ifndef UTOPIA_PLUGIN_FUNCTION_IMPL_HPP
#define UTOPIA_PLUGIN_FUNCTION_IMPL_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

extern "C" {
#include "utopia_plugin_Function.h"
}

namespace utopia {

    class Communicator;

    class PluginFunctionImpl : public Configurable {
    private:
        typedef int (*utopia_plugin_Function_create_crs_graph_t)(const plugin_Function_t *,
                                                                 ptrdiff_t *,
                                                                 ptrdiff_t *,
                                                                 ptrdiff_t *,
                                                                 plugin_idx_t **,
                                                                 plugin_idx_t **);

        typedef int (*utopia_plugin_Function_create_vector_t)(const plugin_Function_t *,
                                                              ptrdiff_t *,
                                                              ptrdiff_t *,
                                                              plugin_scalar_t **);
        typedef int (*utopia_plugin_Function_destroy_vector_t)(const plugin_Function_t *, plugin_scalar_t *values);

        typedef int (*utopia_plugin_Function_value_t)(const plugin_Function_t *,
                                                      const plugin_scalar_t *,
                                                      plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_gradient_t)(const plugin_Function_t *,
                                                         const plugin_scalar_t *const,
                                                         plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_hessian_crs_t)(const plugin_Function_t *,
                                                            const plugin_scalar_t *const,
                                                            const plugin_idx_t *const,
                                                            const plugin_idx_t *const,
                                                            plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_apply_t)(const plugin_Function_t *,
                                                      const plugin_scalar_t *const,
                                                      plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_apply_constraints_t)(const plugin_Function_t *, plugin_scalar_t *const);
        typedef int (*utopia_plugin_Function_apply_zero_constraints_t)(const plugin_Function_t *,
                                                                       plugin_scalar_t *const);
        typedef int (*utopia_plugin_Function_copy_constrained_dofs_t)(const plugin_Function_t *,
                                                                      const plugin_scalar_t *const,
                                                                      plugin_scalar_t *const);

        utopia_plugin_Function_create_crs_graph_t create_crs_graph_{nullptr};

        utopia_plugin_Function_create_vector_t create_vector_{nullptr};
        utopia_plugin_Function_destroy_vector_t destroy_vector_{nullptr};

        utopia_plugin_Function_value_t value_{nullptr};
        utopia_plugin_Function_gradient_t gradient_{nullptr};
        utopia_plugin_Function_hessian_crs_t hessian_crs_{nullptr};
        utopia_plugin_Function_apply_t apply_{nullptr};
        utopia_plugin_Function_apply_constraints_t apply_constraints_{nullptr};
        utopia_plugin_Function_apply_zero_constraints_t apply_zero_constraints_{nullptr};
        utopia_plugin_Function_copy_constrained_dofs_t copy_constrained_dofs_{nullptr};

    public:
        inline int create_crs_graph(ptrdiff_t *nlocal,
                                    ptrdiff_t *nglobal,
                                    ptrdiff_t *nnz,
                                    plugin_idx_t **rowptr,
                                    plugin_idx_t **colidx) const {
            return create_crs_graph_(&info, nlocal, nglobal, nnz, rowptr, colidx);
        }

        inline int create_vector(ptrdiff_t *nlocal, ptrdiff_t *nglobal, plugin_scalar_t **values) const {
            return create_vector_(&info, nlocal, nglobal, values);
        }

        inline int destroy_vector(plugin_scalar_t *values) const { return destroy_vector_(&info, values); }

        inline int value(const plugin_scalar_t *x, plugin_scalar_t *const value) const {
            return value_(&info, x, value);
        }

        inline int gradient(const plugin_scalar_t *const x, plugin_scalar_t *const g) const {
            return gradient_(&info, x, g);
        }

        inline int hessian_crs(const plugin_scalar_t *const x,
                               const plugin_idx_t *const rowptr,
                               const plugin_idx_t *const colidx,
                               plugin_scalar_t *const values) const {
            return hessian_crs_(&info, x, rowptr, colidx, values);
        }

        inline int apply(const plugin_scalar_t *const x, plugin_scalar_t *const y) const { return apply_(&info, x, y); }

        inline int apply_constraints(plugin_scalar_t *const x) const { return apply_constraints_(&info, x); }

        inline int apply_zero_constraints(plugin_scalar_t *const x) const { return apply_zero_constraints_(&info, x); }

        inline int copy_constrained_dofs(const plugin_scalar_t *const src, plugin_scalar_t *const dest) const {
            return copy_constrained_dofs_(&info, src, dest);
        }

        void read(Input &in) override;

        void initialize(const Communicator &comm);

        virtual ~PluginFunctionImpl();

    private:
        plugin_Function_t info;

        typedef int (*utopia_plugin_Function_init_t)(plugin_Function_t *);
        typedef int (*utopia_plugin_Function_destroy_t)(plugin_Function_t *);

        utopia_plugin_Function_init_t init_{nullptr};
        utopia_plugin_Function_destroy_t destroy_{nullptr};

        void *handle{nullptr};
    };

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class PluginFunction : public PluginFunctionImpl {};

}  // namespace utopia

#endif  // UTOPIA_PLUGIN_FUNCTION_IMPL_HPP
