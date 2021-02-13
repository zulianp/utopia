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
            : hessian_approx_strategy_(hessian_approx) {}

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

            std::cout << "---- init again???? \n";
            // NOTE, we are using dummy variables....
            hessian_approx_strategy_->initialize(r, q);
        }

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
                hessian_approx_strategy_->update(alpha_p, q, x, x);

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

                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region4");
                x += alpha * p;
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region4.1");
                r_new = r - alpha * Ap;
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region4.2");
                r_norm = norm2(r_new);
                UTOPIA_NO_ALLOC_END();

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

        bool reset_initial_guess_{false};
        bool initialized_{false};
        bool apply_gradient_descent_step_{false};
        Layout layout_;

        // This fields are not to be copied anywhere
        Vector r, p, q, Ap, alpha_p, r_new, z, z_new;
        std::shared_ptr<HessianApproximation> hessian_approx_strategy_;
        bool hessian_approx_init_{false};
    };
}  // namespace utopia

#endif  // UTOPIA_CONJUGATE_GRADIENT_BFGS_PRECOND_HPP
