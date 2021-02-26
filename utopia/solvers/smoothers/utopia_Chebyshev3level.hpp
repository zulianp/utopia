#ifndef UTOPIA_CHEBYSHEV_HPP
#define UTOPIA_CHEBYSHEV_HPP

#include <memory>
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_PowerMethod.hpp"

namespace utopia {

    /**
     * @brief      Conjugate Gradient solver. Works with all utopia tensor types.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Chebyshev3level final : public OperatorBasedLinearSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        using Preconditioner = utopia::Preconditioner<Vector>;

    public:
        using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;
        using Super::solve;
        using Super::update;

        Chebyshev3level() {}

        Chebyshev3level(const Chebyshev3level &other) : PreconditionedSolverInterface<Vector>(other) {}

        void read(Input &in) override {
            OperatorBasedLinearSolver<Matrix, Vector>::read(in);

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

        bool check_norms() const { return check_norms_; }
        void check_norms(const bool &check_norms) { check_norms_ = check_norms; }

        void init_memory(const Layout &layout) override {
            assert(layout.local_size() > 0);
            OperatorBasedLinearSolver<Matrix, Vector>::init_memory(layout);
            power_method_.init_memory(layout);

            // resets all buffers in case the size has changed
            r_.zeros(layout);
            help_.zeros(layout);
            p_.zeros(layout);

            initialized_ = true;
            layout_ = layout;
        }

        void print_usage(std::ostream &os) const override {
            OperatorBasedLinearSolver<Matrix, Vector>::print_usage(os);
        }

        void update(const Operator<Vector> &A) override {
            const auto layout_rhs = row_layout(A);

            if (!initialized_ || !layout_rhs.same(layout_)) {
                init_memory(layout_rhs);
            }
        }

        Chebyshev3level *clone() const override { return new Chebyshev3level(*this); }

        void copy(const Chebyshev3level &other) { Super::operator=(other); }

        bool smooth(const Vector &rhs, Vector &x) override {
            SizeType temp = this->max_it();
            this->max_it(this->sweeps());
            this->check_norms(false);
            solve(operator_cast<Vector>(*this->get_operator()), rhs, x);
            this->max_it(temp);
            return true;
        }

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override {
            Scalar eigMax = scale_max_eig_ * power_method_.get_max_eig(A, b);
            Scalar eigMin = scale_min_eig_ * eigMax;

            Scalar avg_eig = (eigMax + eigMin) / 2.0;
            Scalar diff_eig = (eigMax - eigMin) / 2.0;

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
                A.apply(p_, help_);
                r_ = r_ + (alpha * help_);

                if (this->check_norms()) {
                    r_norm = norm2(r_);
                } else {
                    r_norm = 9e9;
                }

                converged = this->check_convergence(it, r_norm, 1, 1);

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {r_norm});
                }

                it++;
            }

            return true;
        }

    private:
        bool initialized_{false};
        bool check_norms_{true};
        Layout layout_;
        Scalar scale_max_eig_{1.2}, scale_min_eig_{0.06};

        // This fields are not to be copied anywhere
        Vector r_, p_, help_;
        PowerMethod<Vector> power_method_;
    };
}  // namespace utopia

#endif
