#ifndef UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
#define UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP

#include <utility>

#include "utopia_ConjugateGradient.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_InexactNewtonInterface.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_PrintInfo.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class NewtonInterface : public NonLinearSolver<Vector> {
    public:
        virtual ~NewtonInterface() {}
        virtual bool solve(Function<Matrix, Vector> &fun, Vector &x) = 0;
    };

    template <class Matrix, class Vector>
    class NewtonBase : public NewtonInterface<Matrix, Vector>, public InexactNewtonInterface<Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using Solver = utopia::LinearSolver<Matrix, Vector>;
        using DiffController = utopia::DiffController<Matrix, Vector>;

        using NewtonInterface<Matrix, Vector>::solve;

        NewtonBase(std::shared_ptr<Solver> linear_solver)
            : NewtonInterface<Matrix, Vector>(),
              InexactNewtonInterface<Vector>(),
              linear_solver_(std::move(linear_solver)),
              check_diff_(false) {}

        ~NewtonBase() override = default;

        virtual bool solve(Function_rhs<Matrix, Vector> &fun, Vector &x, const Vector &rhs) {
            fun.set_rhs(rhs);
            bool converged = this->solve(fun, x);
            fun.reset_rhs();
            return converged;
        }

        /**
         * @brief      Enables the differentiation control.
         *
         * @param[in]  checkDiff  Option, if eanable diff_control or no.
         */
        void enable_differentiation_control(bool checkDiff) { check_diff_ = checkDiff; }

        inline bool differentiation_control_enabled() const { return check_diff_; }

        bool check_values(const SizeType iterations,
                          const Function<Matrix, Vector> &fun,
                          const Vector &x,
                          const Vector &gradient,
                          const Matrix &hessian) {
            if (check_diff_ && !controller_.check(fun, x, gradient, hessian)) {
                this->exit_solver(iterations, ConvergenceReason::DIVERGED_INNER);
                return false;
            }

            return true;
        }

        void read(Input &in) override {
            NewtonInterface<Matrix, Vector>::read(in);
            in.get("check_diff", check_diff_);
            in.get_deprecated("check-diff", "check_diff", check_diff_);

            if (check_diff_) {
                in.get("diff_controller", controller_);
                in.get_deprecated("diff-controller", "diff_controller", controller_);
            }

            if (linear_solver_) {
                in.get("linear_solver", *linear_solver_);
                in.get_deprecated("linear-solver", "linear_solver", *linear_solver_);
            }
        }

        void print_usage(std::ostream &os) const override {
            NewtonInterface<Matrix, Vector>::print_usage(os);
            this->print_param_usage(
                os, "check_diff", "bool", "Enables finite difference controller", std::to_string(check_diff_));

            if (linear_solver_) {
                this->print_param_usage(
                    os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "-");
                linear_solver_->print_usage(os);
            } else {
                this->print_param_usage(
                    os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "- (null)");
            }
        }

        /**
         * @brief      Changes linear solver used inside of nonlinear-solver.
         *
         * @param[in]  linear_solver  The linear solver
         */
        virtual void set_linear_solver(const std::shared_ptr<Solver> &linear_solver) { linear_solver_ = linear_solver; }

        inline DiffController &controller() { return controller_; }

    public:
        inline std::shared_ptr<Solver> linear_solver() const { return linear_solver_; }

        void init_memory(const Layout &layout) { linear_solver_->init_memory(layout); }

    protected:
        inline bool linear_solve(const Matrix &mat, const Vector &rhs, Vector &sol) {
            linear_solver_update(mat);
            return linear_solver_apply(rhs, sol);
        }

        inline bool linear_solve(const Matrix &mat, const Matrix &prec, const Vector &rhs, Vector &sol) {
            linear_solver_update(mat, prec);
            return linear_solver_apply(rhs, sol);
        }

        void linear_solver_update(const Matrix &mat) { linear_solver_->update(make_ref(mat)); }
        void linear_solver_update(const Matrix &mat, const Matrix &prec) {
            static_cast<PreconditionedSolver<Matrix, Vector> *>(linear_solver_.get())
                ->update(make_ref(mat), make_ref(prec));
        }

        inline bool linear_solver_apply(const Vector &rhs, Vector &sol) {
            this->solution_status_.num_linear_solves++;
            auto flg = linear_solver_->apply(rhs, sol);
            if (auto *it_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(linear_solver_.get())) {
                auto sol_status_ls = it_solver->solution_status();
                this->solution_status_.sum_linear_its += sol_status_ls.iterates;
            }

            return flg;
        }

        inline bool has_preconditioned_solver() {
            return dynamic_cast<PreconditionedSolver<Matrix, Vector> *>(linear_solver_.get());
        }

        std::shared_ptr<Solver> linear_solver_; /*!< Linear solver parameters. */
        DiffController controller_;
        bool check_diff_; /*!< Enable differentiation control. */
    };
}  // namespace utopia

#endif  // UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
