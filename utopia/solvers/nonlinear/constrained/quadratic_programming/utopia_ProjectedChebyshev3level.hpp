#ifndef UTOPIA_PROJECTED_CHEBYSHEV_HPP
#define UTOPIA_PROJECTED_CHEBYSHEV_HPP

#include <memory>

#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_PowerMethod.hpp"
#include "utopia_QPSolver.hpp"

namespace utopia {

    /**
     * @brief      Conjugate Gradient solver. Works with all utopia tensor types.
     * @tparam     Matrix
     * @tparam     Vector
     */

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class ProjectedChebyshev3level final : public OperatorBasedQPSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using Super = utopia::OperatorBasedQPSolver<Matrix, Vector>;

    public:
        using Super::solve;
        using Super::update;

        ProjectedChebyshev3level() {}

        ProjectedChebyshev3level(const ProjectedChebyshev3level &other)
            : VariableBoundSolverInterface<Vector>(other), PreconditionedSolverInterface<Vector>(other),  Super(other),
              scale_max_eig_(other.scale_max_eig_),
              scale_min_eig_(other.scale_min_eig_),
              power_method_(other.power_method_) {}

        void read(Input &in) override {
            Super::read(in);

            in.get("scale_max_eig", scale_max_eig_);
            in.get("scale_min_eig", scale_min_eig_);

            Scalar eps_eig_est_;
            SizeType power_method_max_it_;
            bool use_rand_vec_init_;

            in.get("power_method_tol", eps_eig_est_);
            in.get("power_method_max_it", power_method_max_it_);
            in.get("use_rand_vec_init", use_rand_vec_init_);

            power_method_.use_rand_vec_init(use_rand_vec_init_);
            power_method_.max_it(power_method_max_it_);
            power_method_.tol(eps_eig_est_);
        }

        void power_method_tol(const Scalar &tol) { power_method_.tol(tol); }
        Scalar power_method_tol() const { return power_method_.tol(); }

        void power_method_max_it(const SizeType &power_method_max_it) { power_method_.max_it(power_method_max_it); }
        SizeType power_method_max_it() const { return power_method_.max_it(); }

        void use_rand_vec_init(const bool &use_rand_vec_init) { power_method_.use_rand_vec_init(use_rand_vec_init); }
        bool use_rand_vec_init() const { return power_method_.use_rand_vec_init(); }

        void scale_max_eig(const Scalar scale_max_eig) { scale_max_eig_ = scale_max_eig; }
        Scalar scale_max_eig() const { return scale_max_eig_; }

        void scale_min_eig(const Scalar scale_min_eig) { scale_min_eig_ = scale_min_eig; }
        Scalar scale_min_eig() const { return scale_min_eig_; }

        void init_memory(const Layout &layout) override {
            assert(layout.local_size() > 0);
            Super::init_memory(layout);
            power_method_.init_memory(layout);

            // resets all buffers in case the size has changed
            r_.zeros(layout);
            help_.zeros(layout);
            p_.zeros(layout);

            initialized_ = true;
            layout_ = layout;
        }

        void print_usage(std::ostream &os) const override { OperatorBasedQPSolver<Matrix, Vector>::print_usage(os); }

        void update(const Operator<Vector> &A) override {
            const auto layout_rhs = row_layout(A);

            if (!initialized_ || !layout_rhs.same(layout_)) {
                init_memory(layout_rhs);
            }

            init_eigs(A);
        }

        ProjectedChebyshev3level *clone() const override { return new ProjectedChebyshev3level(*this); }

        void copy(const ProjectedChebyshev3level &other) {
            Super::operator=(other);
            power_method_ = other.power_method_;
            scale_max_eig_ = other.scale_max_eig_;
            scale_min_eig_ = other.scale_min_eig_;
        }

        bool smooth(const Vector &rhs, Vector &x) override {
            SizeType temp = this->max_it();
            this->max_it(this->sweeps());
            this->norm_frequency(0);
            solve(operator_cast<Vector>(*this->get_operator()), rhs, x);
            this->max_it(temp);
            return true;
        }

        void init_eigs(const Operator<Vector> &A) {
            this->eigMax_ = scale_max_eig_ * power_method_.get_max_eig(A);
            this->eigMin_ = scale_min_eig_ * this->eigMax_;
        }

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override {
            Scalar avg_eig = (this->eigMax_ + this->eigMin_) / 2.0;
            Scalar diff_eig = (this->eigMax_ - this->eigMin_) / 2.0;

            A.apply(x, help_);
            r_ = help_ - b;

            SizeType it = 0;
            Scalar r_norm = 9e9, alpha = 0.0, beta = 0.0;
            bool converged = false;

            while (!converged) {
                if (it == 0) {
                    p_ = -1.0 * r_;
                    alpha = 1.0 / avg_eig;
                } else if (it == 1) {
                    Scalar diff_eig_alpha = (diff_eig * alpha);
                    beta = 0.5 * diff_eig_alpha * diff_eig_alpha;
                    alpha = 1.0 / (avg_eig - (beta / alpha));
                    p_ = -1.0 * r_ + (beta * p_);
                } else {
                    Scalar diff_eig_alpha = (diff_eig * alpha);
                    beta = 0.25 * diff_eig_alpha * diff_eig_alpha;
                    alpha = 1.0 / (avg_eig - (beta / alpha));
                    p_ = -1.0 * r_ + (beta * p_);
                }

                x = x + (alpha * p_);
                this->make_iterate_feasible(x);

                it++;
                converged = this->check_convergence(it, r_norm, 1, 1);
                if (converged) {
                    return true;
                }

                A.apply(x, help_);
                r_ = help_ - b;

                if (this->compute_norm(it) || this->verbose()) {
                    r_norm = norm2(r_);
                } else {
                    r_norm = 9e9;
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {r_norm});
                }
            }

            return true;
        }

    private:
        bool initialized_{false};
        Layout layout_;
        Scalar scale_max_eig_{1.2}, scale_min_eig_{0.06};
        Scalar eigMax_, eigMin_;

        // This fields are not to be copied anywhere
        Vector r_, p_, help_;
        PowerMethod<Vector> power_method_;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTED_CHEBYSHEV_HPP