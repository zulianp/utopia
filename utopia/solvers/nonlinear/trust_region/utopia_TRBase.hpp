#ifndef UTOPIA_SOLVER_TRUSTREGION_BASE_HPP
#define UTOPIA_SOLVER_TRUSTREGION_BASE_HPP

#include <algorithm>
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NumericalTollerance.hpp"
#include "utopia_TRSubproblem.hpp"

namespace utopia {

    template <class Vector>
    class TrustRegionParams : public virtual Configurable {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

    public:
        TrustRegionParams()
            : delta_max_(1e14),
              delta_min_(1e-14),
              delta0_(1.0),
              gamma1_(0.25),
              gamma2_(2.5),
              eta1_(0.05),
              eta2_(0.9),
              rho_tol_(0.005),
              eps_(1e-14) {}

        ~TrustRegionParams() override = default;

        Scalar delta_max() const { return delta_max_; }
        Scalar delta_min() const { return delta_min_; }
        Scalar delta0() const { return delta0_; }
        Scalar gamma1() const { return gamma1_; }
        Scalar gamma2() const { return gamma2_; }
        Scalar eta1() const { return eta1_; }
        Scalar eta2() const { return eta2_; }
        Scalar rho_tol() const { return rho_tol_; }
        Scalar eps() const { return eps_; }

        void delta_max(const Scalar &delta_max_in) { delta_max_ = delta_max_in; };
        void delta_min(const Scalar &delta_min_in) { delta_min_ = delta_min_in; };
        void delta0(const Scalar &delta0_in) { delta0_ = delta0_in; };
        void gamma1(const Scalar &gamma1_in) { gamma1_ = gamma1_in; };
        void gamma2(const Scalar &gamma2_in) { gamma2_ = gamma2_in; };
        void eta1(const Scalar &eta1_in) { eta1_ = eta1_in; };
        void eta2(const Scalar &eta2_in) { eta2_ = eta2_in; };
        void rho_tol(const Scalar &rho_tol_in) { rho_tol_ = rho_tol_in; };
        void eps(const Scalar &eps_in) { eps_ = eps_in; };

        void read(Input &in) override {
            in.get("delta_max", delta_max_);
            in.get("delta_min", delta_min_);
            in.get("delta0", delta0_);
            in.get("gamma1", gamma1_);
            in.get("gamma2", gamma2_);
            in.get("eta1", eta1_);
            in.get("eta2", eta2_);
            in.get("rho_tol", rho_tol_);
            in.get("eps", eps_);
        }

        void print_usage(std::ostream &os) const override {
            this->print_param_usage(os, "delta_max", "real", "Maximum value of tr. radius.", "1e14");
            this->print_param_usage(os, "delta_min", "real", "Minimum value of tr. radius.", "1e-14");
            this->print_param_usage(os, "delta0", "real", "Initial value of tr. radius.", "1.0");

            this->print_param_usage(os, "gamma1", "real", "Factor use to shrink tr. radius.", "0.5");
            this->print_param_usage(os, "gamma2", "real", "Factor use to enlarge tr. radius.", "2.0");

            this->print_param_usage(os, "eta1", "real", "Threshold for rho to shrink tr. radius.", "0.25");
            this->print_param_usage(os, "eta2", "real", "Threshold for rho to enlarge tr. radius.", "0.75");

            this->print_param_usage(os, "rho_tol", "real", "Threshold for rho to take trial point.", "0.005");
            this->print_param_usage(os, "eps", "real", "Numerical tolerance.", "1e-14");
        }

    private:
        Scalar delta_max_;
        Scalar delta_min_;
        Scalar delta0_;
        Scalar gamma1_;
        Scalar gamma2_;
        Scalar eta1_;
        Scalar eta2_;
        Scalar rho_tol_;
        Scalar eps_;
    };

    /**
     * @brief      Base class for all TR solvers. Contains all general routines related to TR solvers.
     *
     */
    template <class Vector>
    class TrustRegionBase : public TrustRegionParams<Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        TrustRegionBase() : TrustRegionParams<Vector>() {}

        ~TrustRegionBase() override = default;

    protected:
        virtual void print_statistics(const SizeType &it, const SizeType &it_successful) {
            std::string path = "log_output_path";
            auto non_data_path = Utopia::instance().get(path);

            if (!non_data_path.empty()) {
                CSVWriter writer{};
                if (mpi_world_rank() == 0) {
                    if (!writer.file_exists(non_data_path)) {
                        writer.open_file(non_data_path);
                        writer.write_table_row<std::string>({"num_its", "it_successful"});
                    } else
                        writer.open_file(non_data_path);

                    writer.write_table_row<Scalar>({Scalar(it), Scalar(it_successful)});
                    writer.close_file();
                }
            }
        }

        virtual Scalar get_pred(const Vector &g, const Operator<Vector> &B, const Vector &p_k) {
            if (empty(Bp_) || size(Bp_) != size(g)) {
                Bp_ = 0.0 * g;
            }

            B.apply(p_k, Bp_);
            return -1.0 * dot(g, p_k) - 0.5 * dot(Bp_, p_k);
        }

        virtual void init_memory(const Layout &layout) { Bp_.zeros(layout); }

        virtual bool check_convergence(Monitor<Vector> &monitor,
                                       const NumericalTollerance<Scalar> &tol,
                                       const SizeType max_it,
                                       const SizeType &it,
                                       const Scalar &g_norm,
                                       const Scalar &r_norm,
                                       const Scalar &s_norm,
                                       const Scalar &delta) const {
            bool converged = false;

            SolutionStatus sol_status;

            // termination because norm of grad is down
            if (g_norm <= tol.absolute_tollerance()) {
                monitor.exit_solver(it, ConvergenceReason::CONVERGED_FNORM_ABS);
                sol_status.reason = ConvergenceReason::CONVERGED_FNORM_ABS;
                converged = true;
            }

            // step size so small that we rather exit than wait for nan's
            else if (s_norm <= tol.step_tollerance()) {
                monitor.exit_solver(it, ConvergenceReason::CONVERGED_SNORM_RELATIVE);
                sol_status.reason = ConvergenceReason::CONVERGED_SNORM_RELATIVE;
                converged = true;
            }

            // step size so small that we rather exit than wait for nan's
            else if (r_norm <= tol.relative_tollerance()) {
                monitor.exit_solver(it, ConvergenceReason::CONVERGED_FNORM_RELATIVE);
                sol_status.reason = ConvergenceReason::CONVERGED_FNORM_RELATIVE;
                converged = true;
            }

            // check number of iterations
            else if (it > max_it) {
                monitor.exit_solver(it, ConvergenceReason::DIVERGED_MAX_IT);
                sol_status.reason = ConvergenceReason::DIVERGED_MAX_IT;
                converged = true;
            }

            // do not hard code this
            else if (delta <= this->delta_min()) {
                monitor.exit_solver(it, ConvergenceReason::CONVERGED_TR_DELTA);
                sol_status.reason = ConvergenceReason::CONVERGED_TR_DELTA;
                converged = true;
            }

            if (converged) {
                sol_status.iterates = it;
                sol_status.gradient_norm = g_norm;
                sol_status.relative_gradient_norm = r_norm;
                sol_status.step_norm = s_norm;

                auto SS_monnitor = monitor.solution_status();
                sol_status.execution_time = SS_monnitor.execution_time;
                monitor.solution_status(sol_status);
            }

            return converged;
        }

        /*!
        \details
                  Trial point acceptance
        @note
        \param rho            - ared/pred
        \param E              - Energy after acceptance
        \param E_k            - current energy
        \param E_k1           - next iterate energy
        \param p_k            - current step
        \param x_k            - curren iterate
        \param x_k1           - new iterate
          */
        virtual bool trial_point_acceptance(const Scalar &rho,
                                            Scalar &E,
                                            const Scalar &E_k,
                                            const Scalar &E_k1,
                                            const Vector &p_k,
                                            const Vector &x_k,
                                            Vector &x_k1) {
            // good reduction, accept trial point
            if (rho >= this->rho_tol()) {
                x_k1 = x_k + p_k;
                E = E_k1;
                return true;
            }
            // otherwise, keep old point
            else {
                x_k1 = x_k;
                E = E_k;
                return false;
            }
        }

        /**
         * @brief      Trial ponit acceptance
         *
         * @param[in]  rho   The rho
         * @param[in]  p_k   The step/direction
         * @param[in]  x_k   The current iterate
         * @param      x_k1  The new_iterate - already in state x_k + p_k
         *
         * @return     accepted or no
         */
        virtual bool trial_point_acceptance(const Scalar &rho, const Vector &x_trial, Vector &x_k) {
            // good reduction, accept trial point
            if (rho >= this->rho_tol()) {
                x_k = x_trial;
                return true;
            }
            // otherwise, keep old point
            else {
                return false;
            }
        }

        /**
         * @brief      Trial ponit acceptance
         *
         * @param[in]  rho   The rho
         *
         * @return     accepted or no
         */
        virtual bool trial_point_acceptance(const Scalar &rho) {
            // good reduction, accept trial point
            if (rho >= this->rho_tol()) {
                return true;
            }
            // otherwise, keep old point
            else {
                return false;
            }
        }

        /*!
        \details
        TR radius update function
        @note
        \param rho            - ared/pred
        \param radius          - tr. radius
        \param p_k            - iterate step
          */
        virtual void delta_update(const Scalar &rho, const Vector &p_k, Scalar &radius, const bool inf_flg = false) {
            if (inf_flg == false) {
                if (rho < this->eta1()) {
                    radius = std::max(Scalar(this->gamma1() * norm2(p_k)), this->delta_min());
                } else if (rho > this->eta2()) {
                    Scalar intermediate = std::max(Scalar(this->gamma2() * norm2(p_k)), radius);
                    radius = std::min(intermediate, this->delta_max());
                }
            } else  // computing update for L_inf norm
            {
                if (rho < this->eta1()) {
                    radius = radius * this->gamma1();
                } else if (rho > this->eta2()) {
                    Scalar intermediate = std::max(Scalar(this->gamma2() * norm_infty(p_k)), radius);

                    // Scalar intermediate = this->gamma2() * radius;
                    radius = std::min(intermediate, this->delta_max());
                }
            }
        }

        virtual void delta_update_inf(const Scalar &rho, const Vector &p_k, Scalar &radius) {
            if (rho >= this->eta2()) {
                radius = std::max(Scalar(this->gamma2() * norm_infty(p_k)), radius);
            } else if (rho < this->eta1() && rho > 0) {
                radius = this->gamma1() * norm_infty(p_k);
            } else if (rho == 0) {
                radius = this->gamma1() * radius;
            }
        }

        // virtual void delta_update_inf(const Scalar &rho, const Vector &p_k, Scalar &radius)
        // {
        //     if(rho >= eta2_)
        //     {
        //       radius = std::max( Scalar(gamma2_ * norm_infty(p_k)), radius);
        //     }
        //     else if (rho < eta1_ )
        //     {
        //       radius = gamma1_ * radius;
        //     }
        // }

        /*!
        \details
        TR radius initialization
        - returns false - choice of tr radius was given by user, alg respect it and doesn't change

        @note
        \param x_k                - initial guess/ current iterate
        \param radius              - tr. radius
          */
        virtual Scalar delta_init(const Vector &x_k, const Scalar &radius0, bool rad_flg) {
            Scalar x_norm = norm2(x_k);
            if (radius0 == 0) {
                rad_flg = true;
                return (x_norm == 0) ? 100 : 100 * x_norm;
            } else {
                rad_flg = false;
                return radius0;
            }
        }

        void read(Input &in) override { TrustRegionParams<Vector>::read(in); }

        void print_usage(std::ostream &os) const override { TrustRegionParams<Vector>::print_usage(os); }

    private:
        Vector Bp_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_TRUSTREGION_BASE_HPP
