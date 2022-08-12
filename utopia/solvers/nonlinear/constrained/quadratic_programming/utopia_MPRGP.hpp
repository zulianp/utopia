#ifndef UTOPIA_CONSTRAINT_QP_MPRGP
#define UTOPIA_CONSTRAINT_QP_MPRGP

#include <string>
#include "utopia_Algorithms.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_QPSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class MPRGP final : public OperatorBasedQPSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using Solver = utopia::LinearSolver<Matrix, Vector>;
        using Super = utopia::OperatorBasedQPSolver<Matrix, Vector>;

    public:
        using Super::solve;
        using Super::update;

        MPRGP() : eps_eig_est_(1e-1), power_method_max_it_(10) {}

        void read(Input &in) override {
            OperatorBasedQPSolver<Matrix, Vector>::read(in);
            in.get("eig_comp_tol", eps_eig_est_);
            in.get("power_method_max_it", power_method_max_it_);
            in.get("hardik_variant", hardik_variant_);
        }

        void print_usage(std::ostream &os) const override {
            OperatorBasedQPSolver<Matrix, Vector>::print_usage(os);
            this->print_param_usage(os, "eig_comp_tol", "double", "Tolerance of eigen solver.", "1e-1");
            this->print_param_usage(
                os, "power_method_max_it", "int", "Maximum number of iterations used inside of power method.", "10");
        }

        MPRGP *clone() const override { return new MPRGP(*this); }

        void update(const Operator<Vector> &A) override {
            UTOPIA_TRACE_REGION_BEGIN("MPRGP::update");

            const auto layout_rhs = row_layout(A);
            if (!initialized_ || !layout_rhs.same(layout_)) {
                init_memory(layout_rhs);
            }

            UTOPIA_TRACE_REGION_END("MPRGP::update");
        }

        bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &x) override {
            UTOPIA_TRACE_REGION_BEGIN("MPRGP::solve(...)");

            if (this->verbose()) {
                this->init_solver("MPRGP comm_size: " + std::to_string(rhs.comm().size()), {"it", "|| g ||"});
            }

            if (this->has_empty_bounds()) {
                this->fill_empty_bounds(layout(rhs));
            } else {
                assert(this->get_box_constraints().valid(layout(rhs)));
            }

            auto &&box = this->get_box_constraints();

            this->update(A);

            // as it is not clear ATM, how to apply preconditioner, we use it at least
            // to obtain initial guess
            if (precond_) {
                // this is unconstrained step
                precond_->apply(rhs, x);
                // projection to feasible set
                this->make_iterate_feasible(x);
            }

            auto ub = box.upper_bound();
            auto lb = box.lower_bound();

            const Scalar gamma = 1.0;
            Scalar alpha_bar = 1;
            if (!hardik_variant_) {
                alpha_bar = 1.95 / this->power_method(A);
            }

            Scalar pAp, beta_beta, fi_fi, gp_dot, g_betta, beta_Abeta;

            SizeType it = 0;
            bool converged = false;
            Scalar gnorm;

            Scalar alpha_cg, alpha_f, beta_sc;

            assert(lb);
            assert(ub);

            this->project(*lb, *ub, x);

            A.apply(x, Ax);
            g = Ax - rhs;

            this->get_fi(x, g, *lb, *ub, fi);
            this->get_beta(x, g, *lb, *ub, beta);

            gp = fi + beta;
            p = fi;

            dots(beta, beta, beta_beta, fi, fi, fi_fi);

            while (!converged) {
                if (beta_beta <= (gamma * gamma * fi_fi)) {
                    A.apply(p, Ap);

                    dots(p, Ap, pAp, g, p, gp_dot);

                    // detecting negative curvature
                    if (pAp <= 0.0) {
                        return true;
                    }

                    alpha_cg = gp_dot / pAp;
                    alpha_f = get_alpha_f(x, p, *lb, *ub, help_f1, help_f2);

                    if (hardik_variant_) {
                        x -= alpha_cg * p;

                        if (alpha_cg <= alpha_f) {
                            g -= alpha_cg * Ap;

                            this->get_fi(x, g, *lb, *ub, fi);
                            beta_sc = dot(fi, Ap) / pAp;
                            p = fi - beta_sc * p;

                        } else {
                            this->project(*lb, *ub, x);
                            A.apply(x, g);
                            g -= rhs;
                            this->get_fi(x, g, *lb, *ub, p);
                        }

                    } else {
                        y = x - alpha_cg * p;

                        if (alpha_cg <= alpha_f) {
                            x = y;
                            g = g - alpha_cg * Ap;
                            this->get_fi(x, g, *lb, *ub, fi);
                            beta_sc = dot(fi, Ap) / pAp;
                            p = fi - beta_sc * p;
                        } else {
                            x = x - alpha_f * p;
                            g = g - alpha_f * Ap;
                            this->get_fi(x, g, *lb, *ub, fi);

                            help_f1 = x - (alpha_bar * fi);
                            this->project(help_f1, *lb, *ub, x);

                            A.apply(x, Ax);
                            g = Ax - rhs;
                            this->get_fi(x, g, *lb, *ub, p);
                        }
                    }
                } else {
                    A.apply(beta, Abeta);

                    dots(g, beta, g_betta, beta, Abeta, beta_Abeta);
                    // detecting negative curvature
                    if (beta_Abeta <= 0.0) {
                        if (this->verbose()) {
                            PrintInfo::print_iter_status(it, {gnorm});
                        }

                        return true;
                    }

                    alpha_cg = g_betta / beta_Abeta;
                    x = x - alpha_cg * beta;
                    g = g - alpha_cg * Abeta;

                    this->get_fi(x, g, *lb, *ub, p);
                }

                this->get_fi(x, g, *lb, *ub, fi);
                this->get_beta(x, g, *lb, *ub, beta);

                gp = fi + beta;

                dots(beta, beta, beta_beta, fi, fi, fi_fi, gp, gp, gnorm);

                gnorm = std::sqrt(gnorm);
                it++;

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {gnorm});
                }

                converged = this->check_convergence(it, gnorm, 1, 1);
            }

            UTOPIA_TRACE_REGION_END("MPRGP::solve(...)");
            return converged;
        }

        void set_eig_comp_tol(const Scalar &eps_eig_est) { eps_eig_est_ = eps_eig_est; }

        void get_fi(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector &fi) const {
            assert(!empty(fi));

            {
                auto d_lb = const_local_view_device(lb);
                auto d_ub = const_local_view_device(ub);
                auto d_x = const_local_view_device(x);
                auto d_g = const_local_view_device(g);
                auto d_fi = local_view_device(fi);

                parallel_for(
                    local_range_device(fi), UTOPIA_LAMBDA(const SizeType i) {
                        // read all
                        const Scalar li = d_lb.get(i);
                        const Scalar ui = d_ub.get(i);
                        const Scalar xi = d_x.get(i);
                        const Scalar gi = d_g.get(i);

                        d_fi.set(i, (li < xi && xi < ui) ? gi : Scalar(0.0));
                    });
            }
        }

        Scalar get_alpha_f(const Vector &x,
                           const Vector &p,
                           const Vector &lb,
                           const Vector &ub,
                           Vector &help_f1,
                           Vector &help_f2) const {
            assert(!empty(help_f1));
            assert(!empty(help_f2));

            {
                auto d_lb = const_local_view_device(lb);
                auto d_ub = const_local_view_device(ub);
                auto d_x = const_local_view_device(x);
                auto d_p = const_local_view_device(p);

                auto h1 = local_view_device(help_f1);
                auto h2 = local_view_device(help_f2);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        // read all for quantities
                        const Scalar li = d_lb.get(i);
                        const Scalar ui = d_ub.get(i);
                        const Scalar xi = d_x.get(i);
                        const Scalar pi = d_p.get(i);

                        // write both helpers
                        h1.set(i, (pi > 0) ? ((xi - li) / pi) : Scalar(1e15));
                        h2.set(i, (pi < 0) ? ((xi - ui) / pi) : Scalar(1e15));
                    });
            }

            return multi_min(help_f1, help_f2);
        }

        void get_beta(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector &beta) const {
            assert(!empty(beta));

            {
                auto d_lb = const_local_view_device(lb);
                auto d_ub = const_local_view_device(ub);
                auto d_x = const_local_view_device(x);
                auto d_g = const_local_view_device(g);
                auto d_beta = local_view_device(beta);

                parallel_for(
                    local_range_device(beta), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar li = d_lb.get(i);
                        const Scalar ui = d_ub.get(i);
                        const Scalar xi = d_x.get(i);
                        const Scalar gi = d_g.get(i);

                        const Scalar val = (device::abs(li - xi) < 1e-14)
                                               ? device::min(0.0, gi)
                                               : ((device::abs(ui - xi) < 1e-14) ? device::max(0.0, gi) : 0.0);

                        d_beta.set(i, val);
                    });
            }
        }

    private:
        Scalar power_method(const Operator<Vector> &A) {
            // Super simple power method to estimate the biggest eigenvalue
            assert(!empty(help_f2));
            help_f2.set(1.0);

            SizeType it = 0;
            bool converged = false;
            Scalar gnorm, lambda = 0.0, lambda_old;

            while (!converged) {
                help_f1 = help_f2;
                A.apply(help_f1, help_f2);
                help_f2 = Scalar(1.0 / Scalar(norm2(help_f2))) * help_f2;

                lambda_old = lambda;

                A.apply(help_f2, help_f1);
                lambda = dot(help_f2, help_f1);

                fi = help_f2 - help_f1;
                gnorm = norm2(fi);

                converged = ((gnorm < eps_eig_est_) || (std::abs(lambda_old - lambda) < eps_eig_est_) ||
                             it > power_method_max_it_)
                                ? true
                                : false;

                it = it + 1;
            }

            if (this->verbose() && mpi_world_rank() == 0)
                utopia::out() << "Power method converged in " << it << " iterations. Largest eig: " << lambda << "  \n";

            return lambda;
        }

    public:
        void init_memory(const Layout &layout) override {
            OperatorBasedQPSolver<Matrix, Vector>::init_memory(layout);

            fi.zeros(layout);
            beta.zeros(layout);
            gp.zeros(layout);
            p.zeros(layout);
            y.zeros(layout);
            Ap.zeros(layout);
            Abeta.zeros(layout);
            Ax.zeros(layout);
            g.zeros(layout);
            help_f1.zeros(layout);
            help_f2.zeros(layout);

            initialized_ = true;
            layout_ = layout;
        }

        void set_preconditioner(const std::shared_ptr<Preconditioner<Vector> > &precond) override {
            precond_ = precond;
        }

    private:
        Vector fi, beta, gp, p, y, Ap, Abeta, Ax, g, help_f1, help_f2;

        Scalar eps_eig_est_;
        SizeType power_method_max_it_;

        bool initialized_{false};
        Layout layout_;

        std::shared_ptr<Preconditioner<Vector> > precond_;
        bool hardik_variant_{false};
    };
}  // namespace utopia

#endif  // UTOPIA_CONSTRAINT_QP_MPRGP
