#ifndef UTOPIA_NL_SOLVE_APP_HPP
#define UTOPIA_NL_SOLVE_APP_HPP

#include "utopia_Main.hpp"
#include "utopia_Multiphysics.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

namespace utopia {

    template <class FunctionSpace>
    class NLSolveApp : public Configurable {
    public:
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using TimeDependentFunction_t = utopia::TimeDependentFunction<FunctionSpace>;
        using Multgrid_t = utopia::SemiGeometricMultigridNew<FunctionSpace>;
        // using FEModelFunction_t = utopia::NewmarkIntegrator<FunctionSpace>;

        // using FEModelFunction_t = ImplicitEulerIntegrator_t;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        // using Newton_t = utopia::Newton<Matrix_t, Vector_t>;
        using Newton_t = utopia::SimpleNewton<Matrix_t, Vector_t>;

        using IO_t = utopia::IO<FunctionSpace>;

        void read(Input &in) override {
            space_ = std::make_shared<FunctionSpace>();
            in.get("space", *space_);

            if (space_->empty()) {
                return;
            }

            std::string integrator;
            in.get("integrator", integrator);

            if (integrator == "Newmark") {
                time_dependent_function_ = utopia::make_unique<NewmarkIntegrator_t>(space_);
                function_ = time_dependent_function_;
            } else if (integrator == "ImplicitEuler") {
                time_dependent_function_ = utopia::make_unique<ImplicitEulerIntegrator_t>(space_);
                function_ = time_dependent_function_;
            } else {
                function_ = std::make_shared<FEModelFunction_t>(space_);
            }

            in.get("problem", *function_);

            bool use_mg = false;
            in.get("use_mg", use_mg);

            if (use_mg) {
                auto mg = std::make_shared<SemiGeometricMultigridNew<FunctionSpace>>();
                mg->set_fine_space(space_);
                linear_solver_ = mg;
            } else {
                linear_solver_ = std::make_shared<OmniLinearSolver_t>();
            }

            solver_ = std::make_shared<Newton_t>(linear_solver_);
            solver_->verbose(true);
            in.get("solver", *solver_);
            in.get("n_time_steps", n_time_steps_);
        }

        bool valid() const { return !space_->empty(); }

        void run_static_problem() {
            Vector_t x;
            function_->create_solution_vector(x);

            if (function_->is_linear()) {
                Matrix_t H;
                Vector_t g;

                function_->hessian_and_gradient(x, H, g);
                g *= -1.0;

                H.convert_to_scalar_matrix();

                linear_solver_->solve(H, g, x);

            } else {
                solver_->solve(*function_, x);
            }

            space_->write("x.e", x);
        }

        void run_time_dependent_problem() {
            Vector_t x;
            time_dependent_function_->create_solution_vector(x);
            time_dependent_function_->setup_IVP(x);

            IO_t io(*space_);
            io.set_output_path("X_t.e");

            for (int t = 0; t < n_time_steps_; ++t) {
                solver_->solve(*time_dependent_function_, x);
                time_dependent_function_->update_IVP(x);

                io.write(x, t, t * time_dependent_function_->delta_time());
                utopia::out() << "time_step: " << t << ", time: " << t * time_dependent_function_->delta_time() << '\n';
            }
        }

        void run() {
            if (function_->is_time_dependent()) {
                run_time_dependent_problem();
            } else {
                run_static_problem();
            }
        }

        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<TimeDependentFunction_t> time_dependent_function_;
        std::shared_ptr<FEFunctionInterface_t> function_;
        std::shared_ptr<LinearSolver_t> linear_solver_;
        std::shared_ptr<Newton_t> solver_;

        int n_time_steps_{2};
    };

}  // namespace utopia

#endif  // UTOPIA_NL_SOLVE_APP_HPP
