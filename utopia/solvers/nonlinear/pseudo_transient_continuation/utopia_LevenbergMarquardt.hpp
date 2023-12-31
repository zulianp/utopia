#ifndef UTOPIA_LEVENBERG_MARQUARDT_HPP
#define UTOPIA_LEVENBERG_MARQUARDT_HPP

#include "utopia_NonlinearLeastSquaresSolver.hpp"

namespace utopia {
    /**
     * @details Implementation of LevenbergMarquardt algorithm chapter 3 of C.T.Kelley, Iterative methods for
     * optimization Damping is performed as described in: An algorithm for least-squares estimation of nonlinear
     * parameters, Marquardt
     */
    template <class Matrix, class Vector>
    class LevenbergMarquardt final : public NonLinearLeastSquaresSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::NonLinearLeastSquaresSolver<Matrix, Vector> NonLinearLeastSquaresSolver;

    public:
        LevenbergMarquardt(const std::shared_ptr<Solver> &linear_solver)
            : NonLinearLeastSquaresSolver(linear_solver),
              omega_up_(2.0),
              omega_down_(0.5),
              mu0_(0.0),
              mu_low_(0.25),
              mu_high_(0.75),
              tau0_(0.001),
              use_diag_scaling_(true) {}

        void omega_up(const Scalar &omega_up) { omega_up_ = omega_up; }
        void omega_down(const Scalar &omega) { omega_down_ = omega; }
        void mu0(const Scalar &mu) { mu0_ = mu; }
        void mu_low(const Scalar &mu) { mu_low_ = mu; }
        void mu_high(const Scalar &mu) { mu_high_ = mu; }
        void tau0(const Scalar &tau) { tau0_ = tau; }
        void use_diag_scaling(const bool use_diag_scaling) { use_diag_scaling_ = use_diag_scaling; }

        Scalar omega_up() const { return omega_up_; }
        Scalar omega_down() const { return omega_down_; }
        Scalar mu0() const { return mu0_; }
        Scalar mu_low() const { return mu_low_; }
        Scalar mu_high() const { return mu_high_; }
        Scalar tau0() const { return tau0_; }
        bool use_diag_scaling() const { return use_diag_scaling_; }

        bool solve(LeastSquaresFunction<Matrix, Vector> &fun, Vector &x_k) override {
            using namespace utopia;
            bool converged = false;

            Scalar it = 0;
            Scalar ared, pred, rho, E_old, E_new, E;

            Vector r_k, p_k, x_trial, g;
            Matrix J_k, J_T, H;

            fun.jacobian(x_k, J_k);
            J_T = transpose(J_k);

            fun.residual(x_k, r_k);
            g = J_T * r_k;
            fun.value(x_k, E_old);

            Matrix I;
            I.identity(square_matrix_layout(layout(r_k)));

            // just to start
            Scalar g0_norm, s_norm, g_norm;
            g0_norm = norm2(r_k);
            g_norm = g0_norm;

            Scalar tau = tau0_;

            if (this->verbose_) {
                this->init_solver("Levenberg-Marquardt", {" it. ", "|| g ||", "J_k", "rho", "tau", "|| p_k ||"});
                PrintInfo::print_iter_status({it, g_norm, E_old, 0.0, tau, 0.0});
            }

            while (!converged) {
                // this line can be done more efficiently
                H = J_T * J_k;

                if (use_diag_scaling_) I = diag(diag(H));

                H += tau * I;

                p_k = 0 * x_k;
                this->linear_solve(H, -1.0 * g, p_k);

                x_trial = x_k + p_k;

                // value of the objective function with correction
                fun.value(x_trial, E_new);
                ared = E_old - E_new;
                pred = get_pred(g, H, p_k);

                rho = ared / pred;
                rho = (std::isfinite(rho)) ? rho : 0.0;
                //----------------------------------------------------------------------------
                //     acceptance of trial point
                //----------------------------------------------------------------------------
                if (rho > mu0_) {
                    x_k = x_trial;

                    fun.residual(x_k, r_k);

                    E_old = E_new;
                    E = E_new;

                    norms2(r_k, p_k, g_norm, s_norm);
                } else {
                    s_norm = norm2(p_k);
                    E = E_old;
                }

                // adjusting damping parameter
                if (rho < mu_low_) {
                    tau = std::max(omega_up_ * tau, tau0_);
                } else if (rho > mu_high_) {
                    tau = omega_down_ * tau;
                }

                tau = (tau < tau0_) ? 0.0 : tau;

                //----------------------------------------------------------------------------
                //    convergence check
                //----------------------------------------------------------------------------
                it++;
                if (this->verbose()) {
                    PrintInfo::print_iter_status({it, g_norm, E, rho, tau, s_norm});
                }

                converged = this->check_convergence(it, g_norm, 9e9, s_norm);

                if (!converged && rho > mu0_) {
                    fun.jacobian(x_k, J_k);
                    J_T = transpose(J_k);

                    g = J_T * r_k;
                }
            }

            return false;
        }

        void read(Input &in) override {
            NonLinearLeastSquaresSolver::read(in);

            in.get("tau0", tau0_);
            in.get("mu_high", mu_high_);
            in.get("mu_low", mu_low_);
            in.get("mu0", mu0_);
            in.get("omega_down", omega_down_);
            in.get("omega_up", omega_up_);
            in.get("use_diag_scaling", use_diag_scaling_);
        }

        void print_usage(std::ostream &os) const override {
            NonLinearLeastSquaresSolver::print_usage(os);

            this->print_param_usage(os, "tau0", "real", "Initial value of pseudo-time-step.", "0.001");

            this->print_param_usage(os, "mu_high", "real", "Treshold for shrinking tau.", "0.75");
            this->print_param_usage(os, "mu_low", "real", "Treshold for enlarging tau.", "0.25");
            this->print_param_usage(os, "mu0", "real", "Treshold for accepting trial point.", "0.0");

            this->print_param_usage(os, "omega_down", "real", "Factor by which we shrink tau.", "0.25");
            this->print_param_usage(os, "omega_up", "real", "Factor by which we enlarge tau.", "2.00");

            this->print_param_usage(
                os, "use_diag_scaling", "bool", "Use diag(J^T J) to scale Hessian. If false, I is used.", "true");
        }

    private:
        Scalar get_pred(const Vector &g, const Matrix &B, const Vector &p_k) {
            return (-1.0 * dot(g, p_k) - 0.5 * dot(B * p_k, p_k));
        }

    private:
        Scalar omega_up_, omega_down_, mu0_, mu_low_, mu_high_, tau0_;
        bool use_diag_scaling_;
    };

}  // namespace utopia

#endif  // UTOPIA_LEVENBERG_MARQUARDT_HPP
