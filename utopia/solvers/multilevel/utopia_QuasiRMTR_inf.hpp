#ifndef UTOPIA_QUASI_RMTR_INF_HPP
#define UTOPIA_QUASI_RMTR_INF_HPP

#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_TRSubproblem.hpp"
#include "utopia_TrustRegionVariableBound.hpp"

#include "utopia_HessianApproximations.hpp"
#include "utopia_RMTR.hpp"

#include "utopia_IdentityTransfer.hpp"
#include "utopia_Transfer.hpp"

namespace utopia {
    /**
     * @brief      The class for Quasi RMTR in infinity norm...
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector, class MLConstraints>
    class QuasiRMTR_inf final : public RMTRBase<Matrix, Vector, FIRST_ORDER_DF>, public MLConstraints {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::RMTRBase<Matrix, Vector, FIRST_ORDER_DF> RMTRBase;
        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

        using TRSubproblem = utopia::MatrixFreeQPSolver<Vector>;
        using TRSubproblemPtr = std::shared_ptr<TRSubproblem>;

        using HessianApproximation = utopia::HessianApproximation<Vector>;
        using HessianApproxPtr = std::shared_ptr<HessianApproximation>;

        typedef utopia::Transfer<Matrix, Vector> Transfer;
        typedef utopia::Level<Matrix, Vector> Level;

        using BoxConstraints = utopia::BoxConstraints<Vector>;
        typedef utopia::RMTRBase<Matrix, Vector, FIRST_ORDER_DF> RMTR;

        using MLConstraints::check_feasibility;

    public:
        QuasiRMTR_inf(const SizeType &n_levels) : RMTR(n_levels), MLConstraints(this->transfer()) {
            hessian_approxs_.resize(n_levels);
        }

        ~QuasiRMTR_inf() override = default;

        void read(Input &in) override {
            RMTRBase::read(in);

            if (tr_subproblems_.size() > 0) {
                in.get("coarse-QPSolver", *tr_subproblems_[0]);

                for (auto i = 1; i < tr_subproblems_.size(); i++) in.get("fine-QPSolver", *tr_subproblems_[i]);
            }

            if (hessian_approxs_.size() > 0) {
                for (auto i = 1; i < hessian_approxs_.size(); i++)
                    in.get("hessian-approx-strategy", *hessian_approxs_[i]);
            }
        }

        void print_usage(std::ostream &os) const override {
            RMTRBase::print_usage(os);

            this->print_param_usage(
                os, "coarse-QPSolver", "MatrixFreeQPSolver", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(
                os, "fine-QPSolver", "MatrixFreeQPSolver", "Input parameters for coarse level QP solver.", "-");
            this->print_param_usage(os,
                                    "hessian-approx-strategy",
                                    "HessianApproximation",
                                    "Input parameters for hessian approximation strategies.",
                                    "-");
        }

        std::string name() override { return "QuasiRMTR_inf"; }

        bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy) {
            if (hessian_approxs_.size() != this->n_levels()) hessian_approxs_.resize(this->n_levels());

            for (auto l = 0; l != hessian_approxs_.size(); ++l)
                hessian_approxs_[l] = std::shared_ptr<HessianApproximation>(strategy->clone());

            return true;
        }

        bool set_hessian_approximation_strategies(const std::vector<HessianApproxPtr> &strategies) {
            if (strategies.size() != this->n_levels()) {
                utopia_error(
                    "utopia::QuasiRMTR::set_hessian_approximation_strategies:: Number of strategies does not equal "
                    "with levels in ML hierarchy. \n");
            }

            hessian_approxs_ = strategies;

            return true;
        }

        bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy,
                                                const SizeType &level) {
            if (hessian_approxs_.size() != this->n_levels()) {
                hessian_approxs_.resize(this->n_levels());
            }

            if (level <= this->n_levels()) {
                hessian_approxs_[level] = strategy;
            } else {
                utopia_error(
                    "utopia::QuasiRMTR::set_tr_strategy:: Requested level exceeds number of levels in ML hierarchy. "
                    "\n");
            }

            return true;
        }

        bool set_coarse_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy) {
            if (tr_subproblems_.size() != this->n_levels()) {
                tr_subproblems_.resize(this->n_levels());
            }

            tr_subproblems_[0] = strategy;
            return true;
        }

        bool set_fine_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy) {
            if (tr_subproblems_.size() != this->n_levels()) {
                tr_subproblems_.resize(this->n_levels());
            }

            // starting from level 1 ....
            for (std::size_t l = 1; l != tr_subproblems_.size(); ++l) {
                tr_subproblems_[l] = std::shared_ptr<TRSubproblem>(strategy->clone());
            }
            return true;
        }

        bool set_tr_strategies(const std::vector<TRSubproblemPtr> &strategies) {
            if (strategies.size() != this->n_levels()) {
                utopia_error(
                    "utopia::RMTR::set_tr_strategies:: Number of tr strategies MUST be equal to number of levels in ML "
                    "hierarchy. \n");
            }

            tr_subproblems_ = strategies;
            return true;
        }

    private:
        void init_memory() override {
            const auto &layouts = this->local_level_layouts();
            bool same_fine_lo = this->init_;

            if (this->init_) {
                same_fine_lo = layouts.back().same(layout(this->memory_.x.back()));
            }

            if (!same_fine_lo) {
                RMTRBase::init_memory();
                MLConstraints::init_memory(layouts);

                for (Scalar l = 0; l < this->n_levels(); l++) {
                    tr_subproblems_[l]->init_memory(layouts[l]);
                    hessian_approxs_[l]->initialize(this->memory_.x[l], this->ml_derivs_.g[l]);
                }
                const SizeType fine_level = this->n_levels() - 1;

                for (auto l = 0; l < fine_level; l++) {
                    this->transfer(l).init_memory();
                }

                this->init_ = true;
            }
        }

        bool check_initialization() override {
            bool flg = RMTRBase::check_initialization();

            if (static_cast<SizeType>(tr_subproblems_.size()) != this->n_levels()) {
                utopia_error("utopia::QuasiRMTR_inf:: number of level QP solvers and levels not equal. \n");
                flg = false;
            }

            if (static_cast<SizeType>(hessian_approxs_.size()) != this->n_levels()) {
                utopia_error("utopia::QuasiRMTR_inf:: number of hessian approxiations and levels do not match. \n");
                flg = false;
            }

            return flg;
        }

        bool delta_update(const Scalar &rho, const SizeType &level, const Vector & /*s_global*/) override {
            Scalar intermediate_delta;

            // we could do also more sophisticated options, but lets not care for the moment ...
            if (rho < this->eta1()) {
                intermediate_delta = std::max(this->gamma1() * this->memory_.delta[level], 1e-15);
            } else if (rho > this->eta2()) {
                intermediate_delta = std::min(this->gamma2() * this->memory_.delta[level], 1e15);
            } else {
                intermediate_delta = this->memory_.delta[level];
            }

            this->memory_.delta[level] = intermediate_delta;

            return false;
        }

        // measuring wrt to feasible set...
        Scalar criticality_measure(const SizeType &level) override {
            return MLConstraints::criticality_measure_inf(level, this->memory_.x[level], this->ml_derivs_.g[level]);
        }

        bool recursion_termination_smoothness(const Vector &g_restricted,
                                              const Vector &g_coarse,
                                              const SizeType &level) override {
            // if we merge calls, reduction can be done together
            Scalar Rg_norm = MLConstraints::criticality_measure_inf(level, this->memory_.x[level], g_restricted);
            Scalar g_norm = MLConstraints::criticality_measure_inf(level, this->memory_.x[level], g_coarse);

            return (Rg_norm >= this->grad_smoothess_termination() * g_norm) ? true : false;
        }

        Scalar get_pred(const SizeType &level) override {
            Scalar l_term = dot(this->ml_derivs_.g[level], this->memory_.s[level]);
            Scalar qp_term = hessian_approxs_[level]->compute_uHu_dot(this->memory_.s[level]);
            return (-l_term - 0.5 * qp_term);
        }

        bool update_level(const SizeType &level) override {
            this->memory_.help[level] = this->ml_derivs_.g[level];
            this->get_multilevel_gradient(this->function(level), level, this->memory_.s_working[level]);
            this->ml_derivs_.y[level] = this->ml_derivs_.g[level] - this->memory_.help[level];

            // swap back....
            this->ml_derivs_.g[level] = this->memory_.help[level];
            hessian_approxs_[level]->update(
                this->memory_.s[level], this->ml_derivs_.y[level], this->memory_.x[level], this->ml_derivs_.g[level]);

            return true;
        }

        void initialize_local_solve(const SizeType &level, const LocalSolveType &solve_type) override {
            // if(!(solve_type == PRE_SMOOTHING && level == this->n_levels()-1))
            // {
            if (solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE) {
                hessian_approxs_[level]->reset();
                // hessian_approxs_[level]->initialize(this->memory_.x[level], this->ml_derivs_.g[level]);
            }
            // }
        }

        bool check_feasibility(const SizeType &level) override {
            return MLConstraints::check_feasibility(level, this->memory_.x[level]);
        }

        void init_level(const SizeType &level) override {
            RMTR::init_level(level);

            const SizeType finer_level = level + 1;
            MLConstraints::init_level(
                level, this->memory_.x[finer_level], this->memory_.x[level], this->memory_.delta[finer_level]);

            // let's see ...
            // this->memory_.delta[level]  = this->delta0();
        }

    public:
        // public because of nvcc
        bool solve_qp_subproblem(const SizeType &level, const bool &flg) override {
            Scalar radius = this->memory_.delta[level];

            // first we need to prepare box of intersection of level constraints with tr. constraints
            std::shared_ptr<Vector> &lb = tr_subproblems_[level]->lower_bound();
            std::shared_ptr<Vector> &ub = tr_subproblems_[level]->upper_bound();

            const Vector &active_lower = this->active_lower(level);
            const Vector &active_upper = this->active_upper(level);

            *lb = active_lower - this->memory_.x[level];
            *ub = active_upper - this->memory_.x[level];

            {
                parallel_transform(*lb, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi)->Scalar {
                    return (xi >= -1.0 * radius) ? xi : -1.0 * radius;
                });

                parallel_transform(*ub, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi)->Scalar {
                    return (xi <= radius) ? xi : radius;
                });
            }

            Scalar atol_level =
                (level == this->n_levels() - 1)
                    ? this->atol()
                    : std::min(this->atol(), this->grad_smoothess_termination() * this->memory_.gnorm[level + 1]);
            if (auto *tr_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(tr_subproblems_[level].get())) {
                if (tr_solver->atol() > atol_level) {
                    tr_solver->atol(atol_level);
                }

                if (flg) {
                    tr_solver->max_it(this->max_QP_coarse_it());
                } else {
                    tr_solver->max_it(this->max_QP_smoothing_it());
                }
            } else {
                assert("QuasiRMTR_inf:: dynamic cas failed. \n");
            }

            this->ml_derivs_.g[level] *= -1.0;
            this->memory_.s[level].set(0.0);
            auto multiplication_action = hessian_approxs_[level]->build_apply_H();
            this->tr_subproblems_[level]->solve(
                *multiplication_action, this->ml_derivs_.g[level], this->memory_.s[level]);
            this->ml_derivs_.g[level] *= -1.0;

            if (has_nan_or_inf(this->memory_.s[level])) {
                this->memory_.s[level].set(0.0);
            } else {
                // ----- just for debugging pourposes, to be commented out in the future...
                MLConstraints::get_projection(*lb, *ub, this->memory_.s[level]);
            }

            return true;
        }

    private:
        std::vector<TRSubproblemPtr> tr_subproblems_;
        std::vector<HessianApproxPtr> hessian_approxs_;
    };

}  // namespace utopia

#endif  // UTOPIA_QUASI_RMTR_INF_HPP