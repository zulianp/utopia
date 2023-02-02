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
        typedef int (*create_crs_graph_t)(const plugin_Function_t *,
                                          ptrdiff_t *,
                                          ptrdiff_t *,
                                          ptrdiff_t *,
                                          plugin_idx_t **,
                                          plugin_idx_t **);

        typedef int (*create_vector_t)(const plugin_Function_t *, ptrdiff_t *, ptrdiff_t *, plugin_scalar_t **);
        typedef int (*destroy_vector_t)(const plugin_Function_t *, plugin_scalar_t *values);

        typedef int (*value_t)(const plugin_Function_t *, const plugin_scalar_t *, plugin_scalar_t *const);

        typedef int (*gradient_t)(const plugin_Function_t *, const plugin_scalar_t *const, plugin_scalar_t *const);

        typedef int (*hessian_crs_t)(const plugin_Function_t *,
                                     const plugin_scalar_t *const,
                                     const plugin_idx_t *const,
                                     const plugin_idx_t *const,
                                     plugin_scalar_t *const);

        typedef int (*apply_t)(const plugin_Function_t *,
                               const plugin_scalar_t *const,
                               const plugin_scalar_t *const,
                               plugin_scalar_t *const);

        typedef int (*apply_constraints_t)(const plugin_Function_t *, plugin_scalar_t *const);
        typedef int (*apply_zero_constraints_t)(const plugin_Function_t *, plugin_scalar_t *const);
        typedef int (*copy_constrained_dofs_t)(const plugin_Function_t *,
                                               const plugin_scalar_t *const,
                                               plugin_scalar_t *const);

        typedef int (*destroy_array_t)(const plugin_Function_t *, void *);
        typedef int (*create_array_t)(const plugin_Function_t *, size_t size, void **);

        typedef int (*report_solution_t)(const plugin_Function_t *, const plugin_scalar_t *const);

        create_crs_graph_t create_crs_graph_{nullptr};

        create_vector_t create_vector_{nullptr};
        destroy_vector_t destroy_vector_{nullptr};

        value_t value_{nullptr};
        gradient_t gradient_{nullptr};
        hessian_crs_t hessian_crs_{nullptr};
        apply_t apply_{nullptr};
        apply_constraints_t apply_constraints_{nullptr};
        apply_zero_constraints_t apply_zero_constraints_{nullptr};
        copy_constrained_dofs_t copy_constrained_dofs_{nullptr};
        destroy_array_t destroy_array_;
        create_array_t create_array_;
        report_solution_t report_solution_;

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

        inline int apply(const plugin_scalar_t *const x,
                         const plugin_scalar_t *const h,
                         plugin_scalar_t *const y) const {
            return apply_(&info, x, h, y);
        }

        inline int apply_constraints(plugin_scalar_t *const x) const { return apply_constraints_(&info, x); }

        inline int apply_zero_constraints(plugin_scalar_t *const x) const { return apply_zero_constraints_(&info, x); }

        inline int copy_constrained_dofs(const plugin_scalar_t *const src, plugin_scalar_t *const dest) const {
            return copy_constrained_dofs_(&info, src, dest);
        }

        void read(Input &in) override;

        void initialize(const Communicator &comm);

        inline int create_array(const size_t size, void **ptr) const { return create_array_(&info, size, ptr); }
        inline int destroy_array(void *ptr) const { return destroy_array_(&info, ptr); }

        inline int report_solution(const plugin_scalar_t *const x) const { return report_solution_(&info, x); }

        virtual ~PluginFunctionImpl();

    private:
        plugin_Function_t info;

        typedef int (*init_t)(plugin_Function_t *);
        typedef int (*destroy_t)(plugin_Function_t *);

        init_t init_{nullptr};
        destroy_t destroy_{nullptr};

        void *handle{nullptr};
    };

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class PluginFunction : public PluginFunctionImpl {};

}  // namespace utopia

#endif  // UTOPIA_PLUGIN_FUNCTION_IMPL_HPP
