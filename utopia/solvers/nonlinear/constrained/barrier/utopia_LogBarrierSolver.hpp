#ifndef UTOPIA_LOG_BARRIER_SOLVER_HPP
#define UTOPIA_LOG_BARRIER_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_Newton.hpp"

#include <iomanip>
#include <limits>

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class LogBarrierSolver final : public NonLinearSolver<Vector>, public VariableBoundSolverInterface<Vector> {
        using NonLinearSolver = utopia::NonLinearSolver<Vector>;
        using VariableBoundSolverInterface = utopia::VariableBoundSolverInterface<Vector>;

        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using Newton = utopia::Newton<Matrix, Vector>;
        using LogBarrierFunction = utopia::LogBarrierFunction<Matrix, Vector>;
        using LSStrategy = utopia::LSStrategy<Vector>;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;

    public:
        LogBarrierSolver(const std::shared_ptr<LinearSolver> &linear_solver =
                             std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE> >())
            : newton_(std::make_shared<Newton>(linear_solver)), function_(std::make_shared<LogBarrierFunction>()) {}

        LogBarrierSolver *clone() const /*override*/ {
            auto ptr = utopia::make_unique<LogBarrierSolver>();
            ptr->newton_ = newton_->clone();
            return ptr.release();
        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) /*override*/ {
            using namespace utopia;
            function_->set_unconstrained_function(make_ref(fun));
            function_->set_box_constraints(make_ref(this->get_box_constraints()));
            function_->reset();

            if (linear_solver_pass_) {
                // One linear solver pass
                Matrix H;
                Vector g;
                fun.hessian_and_gradient(x, H, g);
                newton_->linear_solver()->solve(H, g, x);
            }

            return newton_->solve(*function_, x);
        }

        void read(Input &in) override {
            newton_->read(in);
            function_->read(in);

            Options()
                .add_option("linear_solver_pass",
                            linear_solver_pass_,
                            "Performs a linear solve before integrating the barrier function.")
                .parse(in);
        }

        /**
         * @brief      Sets strategy for computing step-size.
         *
         * @param[in]  strategy  The line-search strategy.
         *
         * @return
         */
        void set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy) {
            newton_->set_line_search_strategy(strategy);
        }

    private:
        std::shared_ptr<Newton> newton_;
        std::shared_ptr<LogBarrierFunction> function_;
        bool linear_solver_pass_{true};
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_SOLVER_HPP
