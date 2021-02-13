#ifndef UTOPIA_CONJUGATE_GRADIENT_BFGS_PRECOND_HPP
#define UTOPIA_CONJUGATE_GRADIENT_BFGS_PRECOND_HPP

#include <memory>
#include "utopia_MatrixFreeLinearSolver.hpp"

namespace utopia {

    /**
     * @brief      Conjugate Gradient solver prejonditioned with Quasi-Newton Hessian approximation of inverse
     *             Implementation follows paper: "Automatic preconditioning by limited memory quasi-Newton updating by
     *             Morales and Nocedal"
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class ConjugateGradientQNPrecond final : public OperatorBasedLinearSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using HessianApproximation = utopia::HessianApproximation<Vector>;

        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        // using Preconditioner = utopia::Preconditioner<Vector>;

    public:
        using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;
        using Super::solve;
        using Super::update;

        ConjugateGradientQNPrecond(const std::shared_ptr<HessianApproximation> &hessian_approx)
            : hessian_approx_strategy_(hessian_approx) {
            if (hessian_approx_strategy_->memory_size() % 2 != 0) {
                utopia_error("Sampling strategy requires m to be devisible by 2");
            }
        }

        void reset_initial_guess(const bool val) { reset_initial_guess_ = val; }
        inline void apply_gradient_descent_step(const bool val) { apply_gradient_descent_step_ = val; }

        void read(Input &in) override {
            Super::read(in);
            in.get("reset_initial_guess", reset_initial_guess_);
            in.get("apply_gradient_descent_step", apply_gradient_descent_step_);
        }

        void init_memory(const Layout &layout) override {
            assert(layout.local_size() > 0);
            OperatorBasedLinearSolver<Matrix, Vector>::init_memory(layout);

            // resets all buffers in case the size has changed
            r.zeros(layout);
            p.zeros(layout);
            q.zeros(layout);
            Ap.zeros(layout);
            alpha_p.zeros(layout);
            r_new.zeros(layout);
            z.zeros(layout);
            z_new.zeros(layout);

            initialized_ = true;
            layout_ = layout;
            hessian_approx_strategy_->initialize(r, q);
            memory_indices_.resize(hessian_approx_strategy_->memory_size());
        }

        void sample_curvature_uniformly(const bool &flg) { uniform_sampling_curvature_ = flg; }
        bool sample_curvature_uniformly() const { return uniform_sampling_curvature_; }

        void print_usage(std::ostream &os) const override {
            OperatorBasedLinearSolver<Matrix, Vector>::print_usage(os);
            this->print_param_usage(
                os, "reset_initial_guess", "bool", "Flag, which decides if initial guess should be reseted.", "false");
        }

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override {
            if (hessian_approx_init_) {
                return preconditioned_solve(A, b, x);
            } else {
                return unpreconditioned_solve(A, b, x);
            }
        }

        void update(const Operator<Vector> &A) override {
            const auto layout_rhs = row_layout(A);

            if (!initialized_ || !layout_rhs.same(layout_)) {
                init_memory(layout_rhs);
            }
        }

        ConjugateGradientQNPrecond *clone() const override { return new ConjugateGradientQNPrecond(*this); }

        void copy(const ConjugateGradientQNPrecond &other) {
            Super::operator=(other);
            reset_initial_guess_ = other.reset_initial_guess_;
        }

        ConjugateGradientQNPrecond(const ConjugateGradientQNPrecond &other)
            : PreconditionedSolverInterface<Vector>(other),
              Super(other),
              reset_initial_guess_(other.reset_initial_guess_),
              apply_gradient_descent_step_(other.apply_gradient_descent_step_) {}

    private:
        bool check_solution(const Operator<Vector> &A, const Vector &x, const Vector &b) const {
            Vector r;
            A.apply(x, r);
            r -= b;

            const Scalar r_norm = norm2(r);

            if (r_norm > 100 * this->atol()) {
                // write("A.m", *this->get_operator());
                // disp(*this->get_operator());
                assert(r_norm <= this->atol());
                return false;
            }

            return true;
        }

        void gradient_descent_step(const Operator<Vector> &A, const Vector &b, Vector &x) {
            A.apply(x, r);
            //-grad = b - A * x
            r = b - r;
            x += r;
        }

        bool unpreconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x) {
            SizeType it = 0;
            Scalar rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9;

            assert(!empty(b));

            // Cheap consistency check
            if (empty(x) || size(x) != size(b)) {
                x.zeros(layout(b));
            } else {
                assert(local_size(x) == local_size(b));
                if (reset_initial_guess_) {
                    x.set(0.);
                }
            }

            if (apply_gradient_descent_step_) {
                gradient_descent_step(A, b, x);
            }

            // r = b - A * x;
            UTOPIA_NO_ALLOC_BEGIN("CG:region1");
            A.apply(x, r);
            r = b - r;
            UTOPIA_NO_ALLOC_END();
            // }

            this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||"});
            bool converged = false;

            SizeType check_norm_each = 1;

            while (!converged) {
                rho = dot(r, r);

                if (rho == 0.) {
                    converged = true;
                    break;
                }

                if (it > 0) {
                    beta = rho / rho_1;
                    UTOPIA_NO_ALLOC_BEGIN("CG:region2");
                    p = r + beta * p;
                    UTOPIA_NO_ALLOC_END();
                } else {
                    UTOPIA_NO_ALLOC_BEGIN("CG:region3");
                    p = r;
                    UTOPIA_NO_ALLOC_END();
                }

                UTOPIA_NO_ALLOC_BEGIN("CG:region4");
                // q = A * p;
                A.apply(p, q);

                Scalar dot_pq = dot(p, q);

                UTOPIA_NO_ALLOC_END();

                if (dot_pq == 0.) {
                    // TODO handle properly
                    utopia_warning("prevented division by zero");
                    converged = true;
                    break;
                }

                UTOPIA_NO_ALLOC_BEGIN("CG:region5");
                alpha = rho / dot_pq;

                // x += alpha * p;
                // r -= alpha * q;

                alpha_p = alpha * p;
                q = alpha * q;
                x += alpha_p;
                r -= q;

                // store hessian approx (x, x is dummy variable here... )
                // hessian_approx_strategy_->update(alpha_p, q, x, x);
                this->update_memory(it, hessian_approx_strategy_->memory_size(), alpha_p, q);

                rho_1 = rho;
                UTOPIA_NO_ALLOC_END();

                if ((it % check_norm_each) == 0) {
                    // r =
                    // A.apply(x, r);
                    // r = b - r;

                    r_norm = norm2(r);

                    if (this->verbose()) {
                        PrintInfo::print_iter_status(it, {r_norm});
                    }

                    converged = this->check_convergence(it, r_norm, 1, 1);
                }

                it++;
            }

            hessian_approx_init_ = true;
            return converged;
        }

        bool preconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x) {
            SizeType it = 0;
            Scalar beta = 0., alpha = 1., r_norm = 9e9;

            z.set(0.0);
            z_new.set(0.0);

            cycle_ = 1;

            // auto precond = this->get_preconditioner();

            if (empty(x) || size(x) != size(b)) {
                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region0");
                x.zeros(layout(b));
            } else {
                assert(local_size(x) == local_size(b));

                if (reset_initial_guess_) {
                    x.set(0.);
                }
            }

            if (apply_gradient_descent_step_) {
                gradient_descent_step(A, b, x);
            }

            UTOPIA_NO_ALLOC_BEGIN("CG_pre:region1");
            A.apply(x, r);
            r = b - r;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("CG_pre:region2");
            // precond->apply(r, z);
            hessian_approx_strategy_->apply_Hinv(r, z);
            p = z;
            UTOPIA_NO_ALLOC_END();

            this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||"});
            bool stop = false;

            while (!stop) {
                // Ap = A*p;
                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region3");
                A.apply(p, Ap);
                alpha = dot(r, z) / dot(p, Ap);
                UTOPIA_NO_ALLOC_END();

                if (std::isinf(alpha) || std::isnan(alpha)) {
                    stop = this->check_convergence(it, r_norm, 1, 1);
                    break;
                }

                alpha_p = alpha * p;
                q = alpha * Ap;

                x += alpha_p;
                r_new = r - q;
                r_norm = norm2(r_new);

                // store hessian approx (x, x is dummy variable here... )
                // hessian_approx_strategy_->update(alpha_p, q, x, x);

                if (r_norm <= this->atol()) {
                    if (this->verbose()) {
                        PrintInfo::print_iter_status(it, {r_norm});
                    }

                    stop = this->check_convergence(it, r_norm, 1, 1);
                    break;
                }

                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region5");
                z_new.set(0.0);
                // precond->apply(r_new, z_new);
                hessian_approx_strategy_->apply_Hinv(r_new, z_new);

                // std::cout << "norm(r_new): " << norm2(r_new) << "  z_new" << norm2(z_new) << "  \n";

                beta = dot(z_new, r_new) / dot(z, r);
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region5.1");
                p = z_new + beta * p;
                r = r_new;
                z = z_new;
                UTOPIA_NO_ALLOC_END();

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {r_norm});
                }

                stop = this->check_convergence(it, r_norm, 1, 1);
                it++;
            }

            if (r_norm <= this->atol()) {
                // FIXME sometimes this fails for some reason
                // assert(check_solution(A, x, b));
                return true;
            } else {
                return false;
            }
        }

        void update_memory(const SizeType &current_iterate, const SizeType &m, const Vector &s, const Vector &y) {
            if (!uniform_sampling_curvature_) {
                hessian_approx_strategy_->update(s, y, s, s);
            } else {
                if (current_iterate < m) {
                    hessian_approx_strategy_->update(s, y, s, s);
                    memory_indices_[current_iterate] = current_iterate;
                } else {
                    Scalar l = current_iterate / (std::pow(2, cycle_)) - m / 2.0 + 1;
                    if ((l == std::floor<SizeType>(l)) && l <= m / 2.0) {
                        SizeType k_prime = (2.0 * l - 1) * std::pow(2.0, cycle_ - 1);
                        SizeType replace_index = findInVector(memory_indices_, k_prime);
                        memory_indices_[replace_index] = current_iterate;
                        // update - which does not exist yet ...
                        hessian_approx_strategy_->replace_at_update_inv(replace_index, s, y);
                        if (l == m / 2.0) {
                            cycle_ += 1;
                        }
                    }
                }
            }
        }

        SizeType findInVector(const std::vector<SizeType> &vector, const SizeType &element) {
            auto it = std::find(vector.begin(), vector.end(), element);
            if (it != vector.end()) {
                return distance(vector.begin(), it);
            } else {
                return -1;
            }
        }

        bool reset_initial_guess_{false};
        bool initialized_{false};
        bool apply_gradient_descent_step_{false};
        Layout layout_;

        // This fields are not to be copied anywhere
        Vector r, p, q, Ap, alpha_p, r_new, z, z_new;
        std::shared_ptr<HessianApproximation> hessian_approx_strategy_;
        bool hessian_approx_init_{false};
        SizeType cycle_{1};

        std::vector<SizeType> memory_indices_;
        bool uniform_sampling_curvature_{true};
    };
}  // namespace utopia

#endif  // UTOPIA_CONJUGATE_GRADIENT_BFGS_PRECOND_HPP
