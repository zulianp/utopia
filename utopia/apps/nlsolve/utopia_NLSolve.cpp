
#include "utopia_AppRunner.hpp"

#include "utopia.hpp"

// #include "utopia_TrivialPreconditioners.hpp"

#ifdef UTOPIA_ENABLE_ISOLVER
#ifdef UTOPIA_ENABLE_PETSC

#include "utopia_plugin_Function_impl.hpp"
using Matrix_t = utopia::PetscMatrix;
using Vector_t = utopia::Traits<Matrix_t>::Vector;
using Scalar_t = utopia::Traits<Vector_t>::Scalar;
using Comm_t = utopia::Traits<Matrix_t>::Communicator;

namespace utopia {

    template <>
    class PluginFunction<Matrix_t> : public Function<Matrix_t, Vector_t>, public Operator<Vector_t> {
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
            isolver_scalar_t *values = 0;

            impl_.create_vector(&nlocal, &nglobal, &values);

            x.wrap(comm_.raw_comm(), nlocal, nglobal, values, [this, values]() { this->impl_.destroy_vector(values); });
        }

        bool apply(const Vector_t &input, Vector_t &output) const override {
            UTOPIA_TRACE_SCOPE("PluginFunction::apply");

            if (!current_solution_) {
                assert(false);
                return false;
            }

            if (output.empty()) {
                output.zeros(layout(input));
            } else {
                output.set(0.);
            }

            {
                auto in_view = local_view_device(input);
                auto out_view = local_view_device(output);
                auto x_view = local_view_device(*current_solution_);

                impl_.apply(&x_view.array()[0], &in_view.array()[0], &out_view.array()[0]);
                impl_.copy_constrained_dofs(&in_view.array()[0], &out_view.array()[0]);
            }

            return true;
        }

        void apply_constraints(Vector_t &x) const {
            auto x_view = local_view_device(x);
            impl_.apply_constraints(&x_view.array()[0]);
        }

        Size size() const override { return {current_solution_->size(), current_solution_->size()}; }
        Size local_size() const override { return {current_solution_->local_size(), current_solution_->local_size()}; }

        // Communicator &comm() override = 0;
        const Communicator &comm() const override { return comm_; }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            UTOPIA_TRACE_SCOPE("PluginFunction::gradient");

            if (g.empty()) {
                g.zeros(layout(x));
            } else {
                g.set(0.);
            }

            {
                auto x_view = local_view_device(x);
                auto g_view = local_view_device(g);

                impl_.gradient(&x_view.array()[0], &g_view.array()[0]);
                impl_.apply_zero_constraints(&g_view.array()[0]);
            }

            double norm_g = norm2(g);
            utopia::out() << "norm_g: " << norm_g << "\n";
            return true;
        }

        bool update(const Vector_t &x) {
            current_solution_ = utopia::make_ref(x);
            return true;
        }

        bool project_onto_feasibile_region(Vector_t &x) const override {
            auto x_view = local_view_device(x);
            impl_.apply_constraints(&x_view.array()[0]);
            return true;
        }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            UTOPIA_TRACE_SCOPE("PluginFunction::hessian");
            // if (H.empty()) {
            ptrdiff_t nlocal;
            ptrdiff_t nglobal;
            ptrdiff_t nnz;
            isolver_idx_t *rowptr;
            isolver_idx_t *colidx;
            isolver_scalar_t *values;

            impl_.create_crs_graph(&nlocal, &nglobal, &nnz, &rowptr, &colidx);
            impl_.create_array(nnz * sizeof(isolver_scalar_t), (void **)&values);
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

            double norm_H = norm2(H);
            utopia::out() << "norm_H: " << norm_H << "\n";

            return true;
        }

        bool hessian(const Vector_t &x, Matrix_t &H, Matrix_t &H_preconditioner) const override { return false; }
        bool report_solution(const Vector_t &x) const {
            auto x_view = local_view_device(x);
            return ISOLVER_FUNCTION_SUCCESS == impl_.report_solution(&x_view.array()[0]);
        }

        Comm_t comm_;
        PluginFunctionImpl impl_;
        std::shared_ptr<const Vector_t> current_solution_;
    };
}  // namespace utopia

void nlsolve(utopia::Input &in) {
    utopia::PluginFunction<Matrix_t> fun;
    fun.read(in);

    Comm_t comm(Comm_t::get_default());
    fun.initialize(comm);

    Vector_t x;
    fun.create_vector(x);
    x.set(0.);

    std::string solver_type = "Newton";
    in.get("solver_type", solver_type);

    if (solver_type == "ConjugateGradient") {
        // Linear
        fun.project_onto_feasibile_region(x);

        Vector_t g;
        fun.create_vector(g);
        fun.gradient(x, g);
        g = -g;

        Vector_t c;
        fun.create_vector(c);
        c.set(0.0);

        utopia::ConjugateGradient<Matrix_t, Vector_t, utopia::HOMEMADE> cg;
        cg.apply_gradient_descent_step(true);
        cg.read(in);

        // x is used for nonlinear material for computing the entries of the Hessian
        fun.update(x);

        bool matrix_free = true;
        in.get("matrix_free", matrix_free);

        if (matrix_free) {
            cg.solve(fun, g, c);
        } else {
            Matrix_t H;
            fun.hessian(x, H);

            // auto prec = std::make_shared<utopia::ILU<Matrix_t, Vector_t>>();
            // utopia::InputParameters params;
            // params.set("block_size", 4);
            // params.set("max_it", 4);
            // prec->read(params);
            // cg.set_preconditioner(prec);

            cg.set_preconditioner(std::make_shared<utopia::InvDiagPreconditioner<Matrix_t, Vector_t>>());

            cg.solve(H, g, c);
        }

        x += c;
    } else if (solver_type == "GradientDescent") {
        auto gd = std::make_shared<utopia::GradientDescent<Vector_t>>();
        gd->damping_parameter(0.9);
        gd->read(in);
        gd->solve(fun, x);
    } else {
        std::shared_ptr<utopia::NewtonBase<Matrix_t, Vector_t>> nlsolver;

        if (solver_type == "Newton") {
            auto newton = std::make_shared<utopia::Newton<Matrix_t>>();
            // auto cg = std::make_shared<utopia::ConjugateGradient<Matrix_t, Vector_t, utopia::HOMEMADE>>();
            // cg->verbose(true);

            // auto ls = std::make_shared<utopia::KSPSolver<Matrix_t, Vector_t>>();
            // ls->pc_type(PCHYPRE);
            // ls->ksp_type(KSPBICG);
            // ls->verbose(true);

            auto ls = std::make_shared<utopia::Factorization<Matrix_t, Vector_t>>();
            newton->set_linear_solver(ls);
            nlsolver = newton;

            bool linear_start = false;
            in.get("linear_start", linear_start);

            if (linear_start) {
                if (!x.comm().rank()) {
                    utopia::out() << "Linear start!\n";
                }
                Matrix_t H;
                Vector_t g;

                fun.hessian(x, H);
                fun.gradient(x, g);
                g = -g;

                fun.apply_constraints(x);
                fun.apply_constraints(g);
                ls->solve(H, g, x);
            } else {
                fun.apply_constraints(x);
            }

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

#endif  // UTOPIA_ENABLE_PETSC
#endif  // UTOPIA_ENABLE_ISOLVER
