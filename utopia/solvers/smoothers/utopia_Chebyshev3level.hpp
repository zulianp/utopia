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

        Chebyshev3level() : eps_eig_est_(1e-5), power_method_max_it_(5) {}

        Chebyshev3level(const Chebyshev3level &other)
            : PreconditionedSolverInterface<Vector>(other),
              eps_eig_est_(other.eps_eig_est_),
              power_method_max_it_(other.power_method_max_it_) {}

        void read(Input &in) override {
            OperatorBasedLinearSolver<Matrix, Vector>::read(in);
            in.get("eig_comp_tol", eps_eig_est_);
            in.get("power_method_max_it", power_method_max_it_);

            in.get("power_method_tol", eps_eig_est_);
            in.get("scale_max_eig", scale_max_eig_);
            in.get("scale_min_eig", scale_min_eig_);
            in.get("use_rand_vec_init", use_rand_vec_init_);
        }

        void power_method_tol(const Scalar &tol) { eps_eig_est_ = tol; }
        Scalar power_method_tol() const { return eps_eig_est_; }

        void power_method_max_it(const SizeType &power_method_max_it) { power_method_max_it_ = power_method_max_it; }
        SizeType power_method_max_it() const { return power_method_max_it_; }

        void use_rand_vec_init(const bool &use_rand_vec_init) { use_rand_vec_init_ = use_rand_vec_init; }
        bool use_rand_vec_init() const { return use_rand_vec_init_; }

        void scale_max_eig(const Scalar scale_max_eig) { scale_max_eig_ = scale_max_eig; }
        Scalar scale_max_eig() const { return scale_max_eig_; }

        void scale_min_eig(const Scalar scale_min_eig) { scale_min_eig_ = scale_min_eig; }
        Scalar scale_min_eig() const { return scale_min_eig_; }

        void init_memory(const Layout &layout) override {
            assert(layout.local_size() > 0);
            OperatorBasedLinearSolver<Matrix, Vector>::init_memory(layout);

            // resets all buffers in case the size has changed
            r_.zeros(layout);
            help_f1.zeros(layout);
            eigenvector_.zeros(layout);
            fi_.zeros(layout);
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

        void copy(const Chebyshev3level &other) {
            Super::operator=(other);
            // reset_initial_guess_ = other.reset_initial_guess_;
        }

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override {
            Scalar eigMax = scale_max_eig_ * this->get_max_eig(A, b);
            Scalar eigMin = scale_min_eig_ * eigMax;

            Scalar avg_eig = (eigMax + eigMin) / 2.0;
            Scalar diff_eig = (eigMax - eigMin) / 2.0;

            A.apply(x, help_f1);
            r_ = help_f1 - b;

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
                    alpha = 1.0 / (avg_eig - beta);
                    p_ = -1.0 * r_ + (beta * p_);
                } else {
                    Scalar diff_eig_alpha = (diff_eig * alpha);
                    beta = 0.25 * diff_eig_alpha;
                    alpha = 1.0 / (avg_eig - beta);
                    p_ = -1.0 * r_ + (beta * p_);
                }

                x = x + (alpha * p_);
                A.apply(p_, help_f1);
                r_ = r_ + (alpha * help_f1);

                r_norm = norm2(r_);

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {r_norm});
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
                it++;
            }

            return true;
        }

    private:
        Scalar get_max_eig(const Operator<Vector> &A, const Vector &rhs) {
            // Super simple power method to estimate the biggest eigenvalue
            assert(!empty(eigenvector_));

            if (use_rand_vec_init_ == true) {
                auto d_eigenvector_ = local_view_device(eigenvector_);

                parallel_for(
                    local_range_device(eigenvector_), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar val = ((Scalar)std::rand() / (RAND_MAX)) + 1;

                        d_eigenvector_.set(i, val);
                    });
            } else if (norm2(eigenvector_) == 0) {
                eigenvector_ = rhs;
            }

            // normalize IG
            eigenvector_ = Scalar(1. / norm2(eigenvector_)) * eigenvector_;

            SizeType it = 0;
            bool converged = false;
            Scalar gnorm, lambda = 0.0, lambda_old;

            while (!converged) {
                help_f1 = eigenvector_;
                A.apply(help_f1, eigenvector_);
                eigenvector_ = Scalar(1.0 / Scalar(norm2(eigenvector_))) * eigenvector_;

                lambda_old = lambda;

                A.apply(eigenvector_, help_f1);
                lambda = dot(eigenvector_, help_f1);

                fi_ = eigenvector_ - help_f1;
                gnorm = norm2(fi_);

                converged = ((gnorm < eps_eig_est_) || (std::abs(lambda_old - lambda) < eps_eig_est_) ||
                             it > power_method_max_it_)
                                ? true
                                : false;

                it = it + 1;
            }

            if (this->verbose())
                utopia::out() << "Power method converged in " << it << " iterations. Largest eig: " << lambda << "  \n";

            return lambda;
        }

        bool initialized_{false};
        Layout layout_;
        Scalar eps_eig_est_;
        SizeType power_method_max_it_;

        Scalar scale_max_eig_{1.2}, scale_min_eig_{0.06};
        bool use_rand_vec_init_{false};

        // This fields are not to be copied anywhere
        Vector r_, p_, help_f1, eigenvector_, fi_;

        PowerMethod<Vector> power_method_;
    };
}  // namespace utopia

#endif
