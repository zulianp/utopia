
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
            return true;
        }
        bool update(const Vector_t &x) { return true; }
        bool project_onto_feasibile_region(Vector_t &x) const override { return true; }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            // if (H.empty()) {
            ptrdiff_t nlocal;
            ptrdiff_t nglobal;
            ptrdiff_t nnz;
            plugin_idx_t *rowptr;
            plugin_idx_t *colidx;
            plugin_scalar_t *values;

            impl_.create_crs_graph(&nlocal, &nglobal, &nnz, &rowptr, &colidx);
            impl_.create_array(nnz * sizeof(plugin_scalar_t), (void **)&values);
            // }

            auto x_view = local_view_device(x);
            impl_.hessian_crs(&x_view.array()[0], rowptr, colidx, values);

            H.wrap(comm_.raw_comm(),
                   nlocal,
                   nlocal,
                   nglobal,
                   nglobal,
                   rowptr,
                   colidx,
                   values,
                   [this, rowptr, colidx, values]() {
                       this->impl_.destroy_array(rowptr);
                       this->impl_.destroy_array(colidx);
                       this->impl_.destroy_array(values);
                   });

            return true;
        }

        bool hessian(const Vector_t &x, Matrix_t &H, Matrix_t &H_preconditioner) const override { return false; }
        bool report_solution(const Vector_t &x) const {
            auto x_view = local_view_device(x);
            return UTOPIA_PLUGIN_SUCCESS == impl_.report_solution(&x_view.array()[0]);
        }

        Comm_t comm_;
        PluginFunctionImpl impl_;
    };
}  // namespace utopia

void nlsolve(utopia::Input &in) {
    utopia::PluginFunction<Matrix_t> fun;
    fun.read(in);

    Comm_t comm(Comm_t::get_default());
    fun.initialize(comm);

    Vector_t x;
    fun.create_vector(x);

    std::string solver_type = "Newton";
    in.get("solver_type", solver_type);

    if (solver_type == "GradientDescent") {
        auto gd = std::make_shared<utopia::GradientDescent<Vector_t>>();
        gd->damping_parameter(0.9);
        gd->read(in);
        gd->solve(fun, x);
    } else {
        std::shared_ptr<utopia::NewtonBase<Matrix_t, Vector_t>> nlsolver;

        if (solver_type == "Newton") {
            auto newton = std::make_shared<utopia::Newton<Matrix_t>>();
            auto cg = std::make_shared<utopia::ConjugateGradient<Matrix_t, Vector_t, utopia::HOMEMADE>>();
            newton->set_linear_solver(cg);
            nlsolver = newton;
        } else {
            auto subproblem = std::make_shared<utopia::SteihaugToint<Matrix_t, Vector_t, utopia::HOMEMADE>>();
            auto trust_region = std::make_shared<utopia::TrustRegion<Matrix_t>>(subproblem);
            nlsolver = trust_region;
        }

        nlsolver->read(in);
        nlsolver->solve(fun, x);
    }

    fun.report_solution(x);
}

// Compile example plug-in
// mpicc -c ../apps/nlsolve/utopia_example_plugin.c -I ../backend/plugin
// mpicc -shared utopia_example_plugin.o -o utopia_example_plugin.dylib
// ./utopia_exec -app nlsolve -path utopia_example_plugin.dylib -damping 0.05
UTOPIA_REGISTER_APP(nlsolve);

#endif
