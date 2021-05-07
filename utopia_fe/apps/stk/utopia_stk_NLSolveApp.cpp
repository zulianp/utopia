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
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        // using Newton_t = utopia::Newton<Matrix_t, Vector_t>;
        using Newton_t = utopia::SimpleNewton<Matrix_t, Vector_t>;

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
        }

        bool valid() const { return !space_->empty(); }

        void run() {
            Vector_t x;
            function_->create_solution_vector(x);
            // This seems to fail with the standard utopia::Newton for some reason
            solver_->solve(*function_, x);
            space_->write("x.e", x);
        }

        std::shared_ptr<FunctionSpace> space_;
        std::unique_ptr<FEModelFunction_t> function_;
        std::shared_ptr<OmniLinearSolver_t> linear_solver_;
        std::shared_ptr<Newton_t> solver_;
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
