#ifndef UTOPIA_PLUGIN_FUNCTION_IMPL_HPP
#define UTOPIA_PLUGIN_FUNCTION_IMPL_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// extern "C" {
#include "isolver_function.h"
// }

namespace utopia {

    class Communicator;

    class PluginFunctionImpl : public Configurable {
    private:
        typedef int (*create_crs_graph_t)(const isolver_function_t *,
                                          ptrdiff_t *,
                                          ptrdiff_t *,
                                          ptrdiff_t *,
                                          isolver_idx_t **,
                                          isolver_idx_t **);

        typedef int (*create_vector_t)(const isolver_function_t *, ptrdiff_t *, ptrdiff_t *, isolver_scalar_t **);
        typedef int (*destroy_vector_t)(const isolver_function_t *, isolver_scalar_t *values);

        typedef int (*value_t)(const isolver_function_t *, const isolver_scalar_t *, isolver_scalar_t *const);

        typedef int (*gradient_t)(const isolver_function_t *, const isolver_scalar_t *const, isolver_scalar_t *const);

        typedef int (*hessian_crs_t)(const isolver_function_t *,
                                     const isolver_scalar_t *const,
                                     const isolver_idx_t *const,
                                     const isolver_idx_t *const,
                                     isolver_scalar_t *const);

        typedef int (*apply_t)(const isolver_function_t *,
                               const isolver_scalar_t *const,
                               const isolver_scalar_t *const,
                               isolver_scalar_t *const);

        typedef int (*apply_constraints_t)(const isolver_function_t *, isolver_scalar_t *const);
        typedef int (*apply_zero_constraints_t)(const isolver_function_t *, isolver_scalar_t *const);
        typedef int (*copy_constrained_dofs_t)(const isolver_function_t *,
                                               const isolver_scalar_t *const,
                                               isolver_scalar_t *const);

        typedef int (*destroy_array_t)(const isolver_function_t *, void *);
        typedef int (*create_array_t)(const isolver_function_t *, size_t size, void **);

        typedef int (*report_solution_t)(const isolver_function_t *, const isolver_scalar_t *const);

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
                                    isolver_idx_t **rowptr,
                                    isolver_idx_t **colidx) const {
            return create_crs_graph_(&info, nlocal, nglobal, nnz, rowptr, colidx);
        }

        inline int create_vector(ptrdiff_t *nlocal, ptrdiff_t *nglobal, isolver_scalar_t **values) const {
            return create_vector_(&info, nlocal, nglobal, values);
        }

        inline int destroy_vector(isolver_scalar_t *values) const { return destroy_vector_(&info, values); }

        inline int value(const isolver_scalar_t *x, isolver_scalar_t *const value) const {
            return value_(&info, x, value);
        }

        inline int gradient(const isolver_scalar_t *const x, isolver_scalar_t *const g) const {
            return gradient_(&info, x, g);
        }

        inline int hessian_crs(const isolver_scalar_t *const x,
                               const isolver_idx_t *const rowptr,
                               const isolver_idx_t *const colidx,
                               isolver_scalar_t *const values) const {
            return hessian_crs_(&info, x, rowptr, colidx, values);
        }

        inline int apply(const isolver_scalar_t *const x,
                         const isolver_scalar_t *const h,
                         isolver_scalar_t *const y) const {
            return apply_(&info, x, h, y);
        }

        inline int apply_constraints(isolver_scalar_t *const x) const { return apply_constraints_(&info, x); }

        inline int apply_zero_constraints(isolver_scalar_t *const x) const { return apply_zero_constraints_(&info, x); }

        inline int copy_constrained_dofs(const isolver_scalar_t *const src, isolver_scalar_t *const dest) const {
            return copy_constrained_dofs_(&info, src, dest);
        }

        void read(Input &in) override;

        void initialize(const Communicator &comm);

        inline int create_array(const size_t size, void **ptr) const { return create_array_(&info, size, ptr); }
        inline int destroy_array(void *ptr) const { return destroy_array_(&info, ptr); }

        inline int report_solution(const isolver_scalar_t *const x) const { return report_solution_(&info, x); }

        virtual ~PluginFunctionImpl();

    private:
        isolver_function_t info;

        typedef int (*init_t)(isolver_function_t *);
        typedef int (*destroy_t)(isolver_function_t *);

        init_t init_{nullptr};
        destroy_t destroy_{nullptr};

        void *handle{nullptr};
    };

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class PluginFunction : public PluginFunctionImpl {};

}  // namespace utopia

#endif  // UTOPIA_PLUGIN_FUNCTION_IMPL_HPP
