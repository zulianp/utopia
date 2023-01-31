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
    public:
        typedef int (*utopia_plugin_Function_create_crs_graph_t)(ptrdiff_t *,
                                                                 ptrdiff_t *,
                                                                 ptrdiff_t *,
                                                                 plugin_idx_t **,
                                                                 plugin_idx_t **);

        typedef int (*utopia_plugin_Function_create_vector_t)(ptrdiff_t *, ptrdiff_t *, plugin_scalar_t **);

        typedef int (*utopia_plugin_Function_value_t)(const plugin_scalar_t *, plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_gradient_t)(const plugin_scalar_t *const, plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_hessian_crs_t)(const plugin_scalar_t *const,
                                                            const plugin_idx_t *const,
                                                            const plugin_idx_t *const,
                                                            plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_apply_t)(const plugin_scalar_t *const, plugin_scalar_t *const);

        typedef int (*utopia_plugin_Function_apply_constraints_t)(const plugin_scalar_t *const);
        typedef int (*utopia_plugin_Function_apply_zero_constraints_t)(const plugin_scalar_t *const);
        typedef int (*utopia_plugin_Function_copy_constrained_dofs_t)(const plugin_scalar_t *const,
                                                                      plugin_scalar_t *const);

        utopia_plugin_Function_create_crs_graph_t create_crs_graph{nullptr};
        utopia_plugin_Function_create_vector_t create_vector{nullptr};
        utopia_plugin_Function_value_t value{nullptr};
        utopia_plugin_Function_gradient_t gradient{nullptr};
        utopia_plugin_Function_hessian_crs_t hessian_crs{nullptr};
        utopia_plugin_Function_apply_t apply{nullptr};
        utopia_plugin_Function_apply_constraints_t apply_constraints{nullptr};
        utopia_plugin_Function_apply_zero_constraints_t apply_zero_constraints{nullptr};
        utopia_plugin_Function_copy_constrained_dofs_t copy_constrained_dofs{nullptr};

        void read(Input &in) override;

        void initialize(Communicator &comm);

        virtual ~PluginFunctionImpl();

    private:
        plugin_Function_t info;

        typedef int (*utopia_plugin_Function_init_t)(MPI_Comm comm, plugin_Function_t *);
        typedef int (*utopia_plugin_Function_destroy_t)(plugin_Function_t *);

        utopia_plugin_Function_init_t init{nullptr};
        utopia_plugin_Function_destroy_t destroy{nullptr};

        void *handle{nullptr};
    };

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class PluginFunction : public PluginFunctionImpl {};

}  // namespace utopia

#endif  // UTOPIA_PLUGIN_FUNCTION_IMPL_HPP
