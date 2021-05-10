#include "utopia_Main.hpp"
#include "utopia_Multiphysics.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

namespace utopia {

    template class NewmarkIntegrator<utopia::stk::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::stk::FunctionSpace>;

    template <class FunctionSpace>
    class NLSolveApp : public Configurable {
    public:
        // using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using FEModelFunction_t = utopia::NewmarkIntegrator<FunctionSpace>;
        // using FEModelFunction_t = ImplicitEulerIntegrator_t;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

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

            function_ = utopia::make_unique<FEModelFunction_t>(space_);
            in.get("problem", *function_);

            linear_solver_ = std::make_shared<OmniLinearSolver_t>();
            solver_ = std::make_shared<Newton_t>(linear_solver_);
            solver_->verbose(true);
            // in.get("solver", *solver_);

            in.get("n_time_steps", n_time_steps_);
        }

        bool valid() const { return !space_->empty(); }

        void run_static_problem() {
            Vector_t x;
            function_->create_solution_vector(x);

            // This seems to fail with the standard utopia::Newton for some reason
            solver_->solve(*function_, x);
            space_->write("x.e", x);
        }

        void run_time_dependent_problem() {
            Vector_t x;
            function_->create_solution_vector(x);
            function_->setup_IVP(x);

            IO_t io(*space_);
            io.set_output_path("X_t.e");

            for (int t = 0; t < n_time_steps_; ++t) {
                solver_->solve(*function_, x);
                function_->update_IVP(x);

                io.write(x, t, t * function_->delta_time());
                utopia::out() << "time_step: " << t << ", time: " << t * function_->delta_time() << '\n';
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
        std::unique_ptr<FEModelFunction_t> function_;
        std::shared_ptr<OmniLinearSolver_t> linear_solver_;
        std::shared_ptr<Newton_t> solver_;

        int n_time_steps_{2};
    };

}  // namespace utopia

void stk_nlsolve(utopia::Input &in) {
    utopia::NLSolveApp<utopia::stk::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "stk_nlsolve: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(stk_nlsolve);
