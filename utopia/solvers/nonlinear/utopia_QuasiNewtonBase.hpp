#ifndef UTOPIA_QUASI_NEWTON_BASE_HPP
#define UTOPIA_QUASI_NEWTON_BASE_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>

namespace utopia {

    template <class Vector>
    class QuasiNewtonBase : public MatrixFreeNonLinearSolver<Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using HessianApproximation = utopia::HessianApproximation<Vector>;
        using MFSolver = utopia::MatrixFreeLinearSolver<Vector>;

        using LSStrategy = utopia::LSStrategy<Vector>;

    public:
        QuasiNewtonBase(const std::shared_ptr<HessianApproximation> &hessian_approx,
                        const std::shared_ptr<MFSolver> &solver)
            : hessian_approx_strategy_(hessian_approx), mf_linear_solver_(solver), alpha_(1.0) {}

        QuasiNewtonBase(const std::shared_ptr<MFSolver> &solver) : mf_linear_solver_(solver), alpha_(1.0) {}

        ~QuasiNewtonBase() override = default;

        void read(Input &in) override {
            MatrixFreeNonLinearSolver<Vector>::read(in);
            in.get("dumping", alpha_);

            if (ls_strategy_) {
                in.get("line-search", *ls_strategy_);
            }
            if (mf_linear_solver_) {
                in.get("linear-solver", *mf_linear_solver_);
            }
            if (hessian_approx_strategy_) {
                in.get("hessian-approx-strategy", *hessian_approx_strategy_);
            }
        }

        void print_usage(std::ostream &os) const override {
            MatrixFreeNonLinearSolver<Vector>::print_usage(os);
            this->print_param_usage(os, "dumping", "double", "Default step size.", "1.0");
            this->print_param_usage(os, "line-search", "LSStrategy", "Input parameters for line-search strategy.", "-");
            this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "-");
            this->print_param_usage(os,
                                    "hessian-approx-strategy",
                                    "HessianApproxStrategy",
                                    "Input parameters for hessian-approximation strategy.",
                                    "-");
        }

        virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy) {
            hessian_approx_strategy_ = strategy;
            return true;
        }

        virtual void set_linear_solver(const std::shared_ptr<MFSolver> &linear_solver) {
            mf_linear_solver_ = linear_solver;
        }

        inline virtual std::shared_ptr<MFSolver> linear_solver() const { return mf_linear_solver_; }

        virtual void update(const Vector &s, const Vector &y, const Vector &x, const Vector &g) {
            if (hessian_approx_strategy_) {
                hessian_approx_strategy_->update(s, y, x, g);
            } else {
                utopia_assert(
                    "QuasiNewtonBase::update, missing hessian approx "
                    "strategy");
            }
        }

        virtual void initialize_approximation(const Vector &x, const Vector &g) {
            if (hessian_approx_strategy_)
                return hessian_approx_strategy_->initialize(x, g);
            else {
                utopia_assert(
                    "QuasiNewtonBase::initialize_approximation, missing hessian approx "
                    "strategy");
            }
        }

        virtual void set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy) { ls_strategy_ = strategy; }

        virtual void dumping_parameter(const Scalar &alpha) { alpha_ = alpha; }

        virtual Scalar dumping_parameter() { return alpha_; }

        virtual void init_memory(const Vector &x, const Vector &g) {
            this->initialize_approximation(x, g);
            auto x_layout = layout(x);

            if (ls_strategy_) {
                ls_strategy_->init_memory(x_layout);
            }

            if (mf_linear_solver_) {
                mf_linear_solver_->init_memory(x_layout);
            }
        }

    protected:
        inline bool linear_solve(const Vector &rhs, Vector &sol) {
            auto multiplication_action = hessian_approx_strategy_->build_apply_H();
            this->solution_status_.num_linear_solves++;
            return mf_linear_solver_->solve(*multiplication_action, rhs, sol);
        }

        inline Scalar get_alpha(FunctionBase<Vector> &fun, const Vector &g, const Vector &x, const Vector &s) {
            if (ls_strategy_) ls_strategy_->get_alpha(fun, g, x, s, alpha_);

            return alpha_;
        }

    protected:
        std::shared_ptr<HessianApproximation> hessian_approx_strategy_;
        std::shared_ptr<MFSolver> mf_linear_solver_;
        Scalar alpha_;                            /*!< Dumping parameter. */
        std::shared_ptr<LSStrategy> ls_strategy_; /*!< Strategy used in order to
                                                     obtain step \f$ \alpha_k \f$ */
    };

}  // namespace utopia
#endif  // UTOPIA_QUASI_NEWTON_BASE_HPP
