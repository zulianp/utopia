
#include "utopia_AppRunner.hpp"

#include "utopia_plugin_Function_impl.hpp"

#include "utopia.hpp"

#ifdef UTOPIA_WITH_PETSC
using Matrix_t = utopia::PetscMatrix;
using Vector_t = utopia::Traits<Matrix_t>::Vector;
using Scalar_t = utopia::Traits<Vector_t>::Scalar;
using Comm_t = utopia::Traits<Matrix_t>::Communicator;

namespace utopia {

    template <>
    class PluginFunction<Matrix_t> : public Function<Matrix_t, Vector_t> {
    public:
        void read(Input &in) override { impl_.read(in); }

        void initialize(const Comm_t &comm) {
            comm_ = comm;
            impl_.initialize(comm_);
        }

        bool value(const Vector_t &x, Scalar &value) const override { return false; }

        void create_vector(Vector_t &x) {
            ptrdiff_t nlocal = 0;
            ptrdiff_t nglobal = 0;
            plugin_scalar_t *values = 0;

            impl_.create_vector(&nlocal, &nglobal, &values);

            x.wrap(comm_.raw_comm(), nlocal, nglobal, values, [this, values]() { this->impl_.destroy_vector(values); });
        }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            if (g.empty()) {
                g.zeros(layout(x));
            }

            auto x_view = local_view_device(x);
            auto g_view = local_view_device(g);

            impl_.gradient(&x_view.array()[0], &g_view.array()[0]);
            return false;
        }
        bool update(const Vector_t &x) { return true; }
        bool project_onto_feasibile_region(Vector_t &x) const override { return true; }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            if (H.empty()) {
                ptrdiff_t nlocal;
                ptrdiff_t nglobal;
                ptrdiff_t nnz;
                plugin_idx_t *rowptr;
                plugin_idx_t *colidx;
                impl_.create_crs_graph(&nlocal, &nglobal, &nnz, &rowptr, &colidx);
            }

            return false;
        }
        bool hessian(const Vector_t &x, Matrix_t &H, Matrix_t &H_preconditioner) const override { return false; }

        Comm_t comm_;
        PluginFunctionImpl impl_;
    };
}  // namespace utopia

void nlsolve(utopia::Input &in) {
    utopia::PluginFunction<Matrix_t> fun;
    fun.read(in);

    Comm_t comm(Comm_t::get_default());
    fun.initialize(comm);

    Vector_t x, g;
    fun.create_vector(x);

    Scalar_t damping = 0.5;
    in.get("damping", damping);

    utopia::GradientDescent<Vector_t> gd;
    gd.dumping_parameter(damping);
    gd.verbose(true);
    gd.solve(fun, x);
    disp(x);
}

// Compile example plug-in
// mpicc -c ../apps/nlsolve/utopia_example_plugin.c -I ../backend/plugin
// mpicc -shared utopia_example_plugin.o -o utopia_example_plugin.dylib
// ./utopia_exec -app nlsolve -path utopia_example_plugin.dylib -damping 0.05
UTOPIA_REGISTER_APP(nlsolve);

#endif
