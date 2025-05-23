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

        Chebyshev3level(const Chebyshev3level &other)
            : PreconditionedSolverInterface<Vector>(other),
              OperatorBasedLinearSolver<Matrix, Vector>(other),
              scale_max_eig_(other.scale_max_eig_),
              scale_min_eig_(other.scale_min_eig_),
              power_method_(other.power_method_) {}

        void read(Input &in) override {
            OperatorBasedLinearSolver<Matrix, Vector>::read(in);

            in.get("scale_max_eig", scale_max_eig_);
            in.get("scale_min_eig", scale_min_eig_);

            Scalar eps_eig_est{1e-2};
            SizeType power_method_max_it{30};
            bool use_rand_vec_init{false};

            in.get("power_method_tol", eps_eig_est);
            in.get("power_method_max_it", power_method_max_it);
            in.get("use_rand_vec_init", use_rand_vec_init);

            power_method_.use_rand_vec_init(use_rand_vec_init);
            power_method_.max_it(power_method_max_it);
            power_method_.tol(eps_eig_est);
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

            init_eigs(A);
        }

        Chebyshev3level *clone() const override { return new Chebyshev3level(*this); }

        void copy(const Chebyshev3level &other) {
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
            UTOPIA_TRACE_SCOPE("Chebyshev3level::solve");
            
            Scalar avg_eig = (this->eigMax_ + this->eigMin_) / 2.0;
            Scalar diff_eig = (this->eigMax_ - this->eigMin_) / 2.0;

            A.apply(x, help_);
            r_ = help_ - b;

            SizeType it = 0;
            Scalar r_norm = 9e9, alpha = 0.0, beta = 0.0;
            bool converged = false;

            this->init_solver("Utopia Chebyshev3level", {"it. ", "||r||"});

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

                it++;
                converged = this->check_convergence(it, r_norm, 1, 1);
                if (converged) {
                    return true;
                }

                A.apply(p_, help_);
                r_ = r_ + (alpha * help_);

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

#endif
