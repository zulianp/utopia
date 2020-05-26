#ifndef UTOPIA_QUASI_TRUST_REGION_HPP
#define UTOPIA_QUASI_TRUST_REGION_HPP

#include "utopia_Dogleg.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_TRBase.hpp"
#include "utopia_TRSubproblem.hpp"

namespace utopia {
    template <class Vector>
    class QuasiTrustRegion final : public TrustRegionBase<Vector>, public QuasiNewtonBase<Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using TRSubproblem = utopia::MatrixFreeTRSubproblem<Vector>;
        using HessianApproximation = utopia::HessianApproximation<Vector>;
        using NonLinearSolver = utopia::QuasiNewtonBase<Vector>;

    public:
        using utopia::QuasiNewtonBase<Vector>::init_memory;

        QuasiTrustRegion(const std::shared_ptr<HessianApproximation> &hessian_approx,
                         const std::shared_ptr<TRSubproblem> &tr_subproblem)
            : NonLinearSolver(hessian_approx, tr_subproblem), initialized_(false) {}

        void read(Input &in) override {
            TrustRegionBase<Vector>::read(in);
            QuasiNewtonBase<Vector>::read(in);
        }

        void print_usage(std::ostream &os) const override {
            TrustRegionBase<Vector>::print_usage(os);
            QuasiNewtonBase<Vector>::print_usage(os);
        }

        /**
         * @brief      QUasi trust region solve.
         *
         * @param      fun   The nonlinear solve function.
         * @param      x_k   Initial gues/ solution
         *
         *
         * @return     true
         */
        bool solve(FunctionBase<Vector> &fun, Vector &x_k) override {
            using namespace utopia;

            // passing solver and parameters into subproblem
            bool converged = false;
            NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

            Scalar delta, product, ared, pred, rho, E_taken, E_old, E_new;  // alpha;

            SizeType it = 0;
            SizeType it_successful = 0;
            Scalar g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity();

            bool rad_flg = false;

            // #define DEBUG_mode

            // TR delta initialization
            delta = this->delta_init(x_k, this->delta0(), rad_flg);

            auto layout_x = layout(x_k);
            if (!initialized_ || !layout_x.same(layout_)) {
                init_memory(layout_x);
            }

            fun.gradient(x_k, g);
            g0_norm = norm2(g);
            g_norm = g0_norm;

            QuasiNewtonBase<Vector>::init_memory(x_k, g);

// print out - just to have idea how we are starting
#ifdef DEBUG_mode
            if (this->verbose_) {
                this->init_solver("QUASI TRUST REGION",
                                  {" it. ",
                                   "|| g ||",
                                   "r_norm",
                                   "<g, dx>",
                                   "J_k",
                                   "J_{k+1/2}",
                                   "J_{k+1}",
                                   "ared",
                                   "pred",
                                   "rho",
                                   "delta_k",
                                   "|| p_k || "});
                PrintInfo::print_iter_status(it, {g_norm});
            }

#else
            if (this->verbose_) {
                this->init_solver("QUASI TRUST REGION",
                                  {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"});
                PrintInfo::print_iter_status(it, {g_norm});
            }
#endif

            it++;

            auto multiplication_action = this->hessian_approx_strategy_->build_apply_H();

            UTOPIA_NO_ALLOC_BEGIN("Quasi TR 1");

            // solve starts here
            while (!converged) {
                fun.value(x_k, E_old);
                //----------------------------------------------------------------------------
                //     new step p_k w.r. ||p_k|| <= delta
                //----------------------------------------------------------------------------
                if (auto *tr_subproblem = dynamic_cast<TRSubproblem *>(this->linear_solver().get())) {
                    p_k.set(0);
                    tr_subproblem->current_radius(delta);
                    g_help = -1.0 * g;
                    tr_subproblem->solve(*multiplication_action, g_help, p_k);
                    this->solution_status_.num_linear_solves++;
                } else {
                    utopia_warning("TrustRegion::Set suitable TR subproblem.... \n ");
                }

                x_trial = x_k + p_k;
                pred = this->get_pred(g, *multiplication_action, p_k);
                //----------------------------------------------------------------------------
                //----------------------------------------------------------------------------
                if (it == 1 && rad_flg) {
                    delta = norm2(p_k);
                    delta *= 0.2;
                }

                // value of the objective function with correction
                fun.value(x_trial, E_new);
                product = dot(g, p_k);  // just to do tests

                // decrease ratio
                ared = E_old - E_new;  // reduction observed on objective function
                pred = std::abs(pred);
                rho = ared / pred;  // decrease ratio

                //----------------------------------------------------------------------------
                //     acceptance of trial point
                //----------------------------------------------------------------------------
                if (ared < 0 || pred < 0) {
                    rho = 0;
                } else if (ared == pred) {
                    rho = 1;
                }

                if (rho >= this->rho_tol()) it_successful++;

                // good reduction, accept trial point
                if (rho >= this->rho_tol()) {
                    x_k += p_k;

                    y = g;
                    fun.gradient(x_k, g);
                    y = g - y;

                    E_taken = E_new;

                }
                // otherwise, keep old point
                else {
                    // Vector grad_trial;
                    fun.gradient(x_trial, g_help);
                    y = g_help - g;

                    E_taken = E_old;
                }

                this->update(p_k, y, x_k, g);

                //----------------------------------------------------------------------------
                //    convergence check
                //----------------------------------------------------------------------------
                norms2(g, p_k, g_norm, s_norm);
                r_norm = g_norm / g0_norm;

#ifdef DEBUG_mode
                if (this->verbose_)
                    PrintInfo::print_iter_status(
                        it, {g_norm, r_norm, product, E_taken, E_old, E_new, ared, pred, rho, delta, s_norm});
#else
                if (this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, E_taken, E_new, rho, delta, s_norm});
#endif

                converged = TrustRegionBase<Vector>::check_convergence(
                    *this, tol, this->max_it(), it, g_norm, r_norm, 9e9, delta);
                //----------------------------------------------------------------------------
                //      tr. radius update
                //----------------------------------------------------------------------------
                this->delta_update(rho, p_k, delta);
                it++;
            }

            UTOPIA_NO_ALLOC_END();

            return true;
        }

        void set_trust_region_strategy(const std::shared_ptr<TRSubproblem> &tr_linear_solver) {
            NonLinearSolver::set_linear_solver(tr_linear_solver);
        }

    private:
        void init_memory(const Layout &layout) override {
            y.zeros(layout);
            p_k.zeros(layout);
            x_trial.zeros(layout);
            g_help.zeros(layout);
            g.zeros(layout);

            TrustRegionBase<Vector>::init_memory(layout);

            initialized_ = true;
            layout_ = layout;
        }

    private:
        Vector g, g_help, y, p_k, x_trial;
        bool initialized_;
        Layout layout_;
    };

}  // namespace utopia

#endif  // UTOPIA_QUASI_TRUST_REGION_HPP
