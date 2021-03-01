#ifndef UTOPIA_CONJUGATE_GRADIENT_BFGS_PRECOND_HPP
#define UTOPIA_CONJUGATE_GRADIENT_BFGS_PRECOND_HPP

#include <memory>
#include "utopia_MatrixFreeLinearSolver.hpp"

namespace utopia {

    /**
     * @brief      Conjugate Gradient solver prejonditioned with Quasi-Newton Hessian approximation of inverse
     *             Implementation follows paper: "Automatic preconditioning by limited memory quasi-Newton updating by
     *             Morales and Nocedal"
     *
     *             Uses heuristic to shift spectrum, if matrix is negative definite
     *
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

    public:
        using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;
        using Super::solve;
        using Super::update;

        ConjugateGradientQNPrecond(const std::shared_ptr<HessianApproximation> &hessian_approx)
            : hessian_approx_strategy_new_(hessian_approx) {
            hessian_approx_strategy_old_ = std::shared_ptr<HessianApproximation>(hessian_approx_strategy_new_->clone());

            if (hessian_approx_strategy_new_->memory_size() % 2 != 0) {
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
            hessian_approx_strategy_new_->initialize(r, q);
            hessian_approx_strategy_old_->initialize(r, q);
            memory_indices_.resize(hessian_approx_strategy_new_->memory_size());
        }

        void sample_curvature_uniformly(const bool &flg) { uniform_sampling_curvature_ = flg; }
        bool sample_curvature_uniformly() const { return uniform_sampling_curvature_; }

        void sample_only_once(const bool &flg) { sample_only_once_ = flg; }
        bool sample_only_once() const { return sample_only_once_; }

        void print_usage(std::ostream &os) const override {
            OperatorBasedLinearSolver<Matrix, Vector>::print_usage(os);
            this->print_param_usage(
                os, "reset_initial_guess", "bool", "Flag, which decides if initial guess should be reseted.", "false");
        }

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override {
            Scalar pAp_pp = 0.0;
            bool flg;

            if (hessian_approx_init_) {
                flg = preconditioned_solve(A, b, x, pAp_pp);
            } else {
                flg = unpreconditioned_solve(A, b, x, pAp_pp);
            }

            while (pAp_pp != 0.0) {
                shift_ = 5.0 * device::max(-pAp_pp, shift_);
                utopia::out() << "------ shifting eig: " << shift_ << " \n";
                flg = preconditioned_solve(A, b, x, pAp_pp);
            }

            return flg;
        }

        void update(const Operator<Vector> &A) override {
            // utopia::out() << "--- update coarse ----- \n";

            const auto layout_rhs = row_layout(A);

            if (!initialized_ || !layout_rhs.same(layout_)) {
                init_memory(layout_rhs);
            }

            shift_ = 0.0;

            // reset precond as operator changed
            if (reset_precond_update_) {
                hessian_approx_init_ = false;
            }
        }

        ConjugateGradientQNPrecond *clone() const override { return new ConjugateGradientQNPrecond(*this); }

        void copy(const ConjugateGradientQNPrecond &other) {
            Super::operator=(other);
            reset_initial_guess_ = other.reset_initial_guess_;
            apply_gradient_descent_step_ = other.apply_gradient_descent_step_;

            // check if only settings get copied, or whole vecs...
            hessian_approx_strategy_new_ = other.hessian_approx_strategy_new_;
            hessian_approx_strategy_old_ = other.hessian_approx_strategy_old_;
            uniform_sampling_curvature_ = other.uniform_sampling_curvature_;
            sample_only_once_ = other.sample_only_once_;
            reset_precond_update_ = other.reset_precond_update_;
        }

        ConjugateGradientQNPrecond(const ConjugateGradientQNPrecond &other)
            : PreconditionedSolverInterface<Vector>(other),
              Super(other),
              reset_initial_guess_(other.reset_initial_guess_),
              apply_gradient_descent_step_(other.apply_gradient_descent_step_),
              hessian_approx_strategy_new_(other.hessian_approx_strategy_new_),
              hessian_approx_strategy_old_(other.hessian_approx_strategy_old_),
              uniform_sampling_curvature_(other.uniform_sampling_curvature_),
              sample_only_once_(other.sample_only_once_),
              reset_precond_update_(other.reset_precond_update_) {}

        bool reset_precond_update() const { return reset_precond_update_; }
        void reset_precond_update(const bool &flg) { reset_precond_update_ = flg; }

    private:
        bool check_solution(const Operator<Vector> &A, const Vector &x, const Vector &b) const {
            Vector r;
            A.apply(x, r);
            if (shift_ != 0.0) {
                r += shift_ * x;
            }

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
            if (shift_ != 0.0) {
                r += shift_ * x;
            }

            //-grad = b - A * x
            r = b - r;
            x += r;
        }

        bool unpreconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x, Scalar &pAp_pp) {
            SizeType it = 0;
            Scalar rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9;

            pAp_pp = 0.0;

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

            A.apply(x, r);
            if (shift_ != 0.0) {
                r += shift_ * x;
            }
            r = b - r;

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
                    p = r + beta * p;
                } else {
                    p = r;
                }

                // q = A * p;
                A.apply(p, q);
                if (shift_ != 0.0) {
                    q += shift_ * p;
                }

                Scalar dot_pq = dot(p, q);

                // checking for negative curvature
                if (dot_pq <= 0.0) {
                    this->check_convergence(it, r_norm, 1, 1);
                    gradient_descent_step(A, b, x);
                    utopia::out() << "---- negative curvature check failed --- \n";
                    pAp_pp = dot_pq / dot(p, p);
                    break;
                }

                alpha = rho / dot_pq;
                alpha_p = alpha * p;
                q = alpha * q;
                x += alpha_p;
                r -= q;

                // store hessian approx (x, x is dummy variable here... )
                // hessian_approx_strategy_->update(alpha_p, q, x, x);
                this->update_memory(it, hessian_approx_strategy_new_->memory_size(), alpha_p, q);

                rho_1 = rho;

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
            if (sample_only_once_) {
                hessian_approx_strategy_old_ =
                    std::shared_ptr<HessianApproximation>(hessian_approx_strategy_new_->clone());
            }

            // std::cout << "max_it_coarse_grid: " << it << "   \n";

            return converged;
        }

        bool preconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x, Scalar &pAp_pp) {
            SizeType it = 0;
            Scalar beta = 0., alpha = 1., r_norm = 9e9;
            Scalar pAp, rz;

            pAp_pp = 0.0;

            z.set(0.0);
            z_new.set(0.0);

            cycle_ = 1;
            if (!sample_only_once_) {
                hessian_approx_strategy_old_ =
                    std::shared_ptr<HessianApproximation>(hessian_approx_strategy_new_->clone());
            }

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

            A.apply(x, r);
            if (shift_ != 0.0) {
                r += shift_ * x;
            }

            r = b - r;

            hessian_approx_strategy_old_->apply_Hinv(r, z);
            p = z;

            this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||"});
            bool stop = false;

            while (!stop) {
                A.apply(p, Ap);
                if (shift_ != 0.0) {
                    Ap += shift_ * p;
                }

                dots(r, z, rz, p, Ap, pAp);
                alpha = rz / pAp;

                if (std::isinf(alpha) || std::isnan(alpha)) {
                    stop = this->check_convergence(it, r_norm, 1, 1);
                    gradient_descent_step(A, b, x);
                    utopia::out() << "---- nan --- \n";
                    break;
                }

                // checking negative curvature
                if (pAp <= 0.0) {
                    stop = this->check_convergence(it, r_norm, 1, 1);
                    gradient_descent_step(A, b, x);
                    utopia::out() << "---- negative curvature check failed --- \n";
                    pAp_pp = pAp / dot(p, p);
                    break;
                }

                alpha_p = alpha * p;
                q = alpha * Ap;

                x += alpha_p;
                r_new = r - q;
                r_norm = norm2(r_new);

                // store hessian approx (x, x is dummy variable here... )
                if (!sample_only_once_) {
                    this->update_memory(it, hessian_approx_strategy_new_->memory_size(), alpha_p, q);
                }

                if (r_norm <= this->atol()) {
                    if (this->verbose()) {
                        PrintInfo::print_iter_status(it, {r_norm});
                    }

                    stop = this->check_convergence(it, r_norm, 1, 1);
                    break;
                }

                z_new.set(0.0);
                hessian_approx_strategy_old_->apply_Hinv(r_new, z_new);

                beta = dot(z_new, r_new) / dot(z, r);

                p = z_new + beta * p;
                r = r_new;
                z = z_new;

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {r_norm});
                }

                stop = this->check_convergence(it, r_norm, 1, 1);

                it++;
            }

            // std::cout << "max_it_coarse_grid: " << it << "   \n";

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
                hessian_approx_strategy_new_->update(s, y, s, s);
            } else {
                if (current_iterate < m) {
                    hessian_approx_strategy_new_->update(s, y, s, s);
                    memory_indices_[current_iterate] = current_iterate;
                } else {
                    Scalar l = current_iterate / (std::pow(2, cycle_)) - m / 2.0 + 1;
                    if ((l == std::floor<SizeType>(l)) && l <= m / 2.0) {
                        SizeType k_prime = (2.0 * l - 1) * std::pow(2.0, cycle_ - 1);
                        SizeType replace_index = findInVector(memory_indices_, k_prime);
                        memory_indices_[replace_index] = current_iterate;
                        // update - which does not exist yet ...
                        hessian_approx_strategy_new_->replace_at_update_inv(replace_index, s, y);
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
        std::shared_ptr<HessianApproximation> hessian_approx_strategy_new_;
        std::shared_ptr<HessianApproximation> hessian_approx_strategy_old_;

        bool hessian_approx_init_{false};
        SizeType cycle_{1};

        std::vector<SizeType> memory_indices_;
        bool uniform_sampling_curvature_{true};
        bool sample_only_once_{true};

        bool reset_precond_update_{true};

        Scalar shift_{0.0};
    };
}  // namespace utopia

#endif  // UTOPIA_CONJUGATE_GRADIENT_BFGS_PRECOND_HPP
