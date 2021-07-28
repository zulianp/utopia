#ifndef UTOPIA_NL_SOLVE_APP_HPP
#define UTOPIA_NL_SOLVE_APP_HPP

#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"
#include "utopia_VelocityImplicitEulerIntegrator.hpp"
#include "utopia_VelocityNewmarkIntegrator.hpp"

#include "utopia_NLSolve.hpp"

namespace utopia {

    template <class FunctionSpace>
    class NLSolveApp : public Configurable {
    public:
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using VelocityImplicitEulerIntegrator_t = utopia::VelocityImplicitEulerIntegrator<FunctionSpace>;

        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using VelocityNewmarkIntegratorIntegrator_t = utopia::VelocityNewmarkIntegrator<FunctionSpace>;

        using TimeDependentFunction_t = utopia::TimeDependentFunction<FunctionSpace>;
        using Multgrid_t = utopia::SemiGeometricMultigridNew<FunctionSpace>;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Newton_t = utopia::Newton<Matrix_t, Vector_t>;

        using IO_t = utopia::IO<FunctionSpace>;

        void read(Input &in) override {
            std::shared_ptr<TimeDependentFunction_t> time_dependent_function;
            std::shared_ptr<FEFunctionInterface_t> function;
            std::shared_ptr<LinearSolver_t> linear_solver;
            std::shared_ptr<Newton_t> solver;

            auto space_ = std::make_shared<FunctionSpace>();
            in.require("space", *space_);

            if (space_->empty()) {
                return;
            }

            std::shared_ptr<FEFunctionInterface<FunctionSpace>> problem;
            std::string functiontype = "";
            in.get("problem", [&functiontype](Input &node) { node.get("type", functiontype); });

            problem = std::make_shared<FEModelFunction_t>(space_);

            std::string integrator;
            in.get("integrator", integrator);

            if (integrator == "Newmark") {
                time_dependent_function = utopia::make_unique<NewmarkIntegrator_t>(problem);
                function = time_dependent_function;
            } else if (integrator == "ImplicitEuler") {
                time_dependent_function = utopia::make_unique<ImplicitEulerIntegrator_t>(problem);
                function = time_dependent_function;
            } else if (integrator == "VelocityImplicitEuler") {
                time_dependent_function = utopia::make_unique<VelocityImplicitEulerIntegrator_t>(problem);
                function = time_dependent_function;
            } else if (integrator == "VelocityNewmark") {
                time_dependent_function = utopia::make_unique<VelocityNewmarkIntegratorIntegrator_t>(problem);
                function = time_dependent_function;
            } else {
                function = problem;
            }

            in.require("problem", *function);

            bool use_mg = false;
            in.get("use_mg", use_mg);

            if (use_mg) {
                auto mg = std::make_shared<SemiGeometricMultigridNew<FunctionSpace>>();
                mg->set_fine_space(space_);
                linear_solver = mg;
            } else {
                linear_solver = std::make_shared<OmniLinearSolver_t>();
            }

            solver = std::make_shared<Newton_t>(linear_solver);
            solver->verbose(true);
            in.get("solver", *solver);

            nlsolve_.init(function);
            nlsolve_.set_solver(solver);
            valid_ = true;
        }

        bool valid() const { return valid_; }
        void run() { nlsolve_.solve(); }

    private:
        NLSolve<FunctionSpace> nlsolve_;
        bool valid_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_NL_SOLVE_APP_HPP
