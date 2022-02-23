#ifndef UTOPIA_GRADIENT_DESCENT_SOLVER_HPP
#define UTOPIA_GRADIENT_DESCENT_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>

namespace utopia {
    /**
     * @brief      The Gradient Descent solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Vector>
    class GradientDescent final : public MatrixFreeNonLinearSolver<Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        using LSStrategy = utopia::LSStrategy<Vector>;

    public:
        GradientDescent() : MatrixFreeNonLinearSolver<Vector>(), alpha_(1.0) {}

        bool solve(FunctionBase<Vector> &fun, Vector &x) override {
            using namespace utopia;

            Vector grad, step;

            Scalar g_norm, g0_norm, r_norm, s_norm;
            SizeType it = 0;

            bool converged = false;

            // notify listener
            fun.update(x);

            fun.gradient(x, grad);
            g0_norm = norm2(grad);
            g_norm = g0_norm;

            this->init_solver("GradientDescent", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha_k"});

            if (this->verbose_) PrintInfo::print_iter_status(it, {g_norm, 1, 0});
            it++;

            while (!converged) {
                // find direction step
                step = -grad;
                if (ls_strategy_) {
                    ls_strategy_->get_alpha(fun, grad, x, step, alpha_);
                    x += alpha_ * step;

                } else {
                    // update x
                    if (fabs(alpha_ - 1) < std::numeric_limits<Scalar>::epsilon()) {
                        x += step;
                    } else {
                        x += alpha_ * step;
                    }
                }

                // notify listener
                fun.update(x);
                fun.gradient(x, grad);

                // norms needed for convergence check
                norms2(grad, step, g_norm, s_norm);
                r_norm = g_norm / g0_norm;

                // // print iteration status on every iteration
                if (this->verbose_) PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm, alpha_});

                // // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, r_norm, s_norm);
                it++;
            }

            this->print_statistics(it);

            return true;
        }

        void read(Input &in) override {
            MatrixFreeNonLinearSolver<Vector>::read(in);
            in.get("dumping", alpha_);

            if (ls_strategy_) {
                in.get("line-search", *ls_strategy_);
            }
        }

        void print_usage(std::ostream &os) const override {
            MatrixFreeNonLinearSolver<Vector>::print_usage(os);
            this->print_param_usage(os, "dumping", "real", "Default step size.", "1.0");
            this->print_param_usage(os, "line-search", "LSStrategy", "Input parameters for line-search strategy.", "-");
        }

        bool set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy) {
            ls_strategy_ = strategy;
            return true;
        }

        void dumping_parameter(const Scalar &alpha) { alpha_ = alpha; }

        const Scalar dumping_parameter() { return alpha_; }

    private:
        Scalar alpha_;                            /*!< Dumping parameter. */
        std::shared_ptr<LSStrategy> ls_strategy_; /*!< Strategy used in order to
                                                     obtain step \f$ \alpha_k \f$ */
    };

}  // namespace utopia
#endif  // UTOPIA_GRADIENT_DESCENT_SOLVER_HPP
