#ifndef UTOPIA_SOLVER_SEMISMOOTH_NEWTON_OLD_HPP
#define UTOPIA_SOLVER_SEMISMOOTH_NEWTON_OLD_HPP

#include "utopia_Core.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_SemismoothNewton.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class Matrix, class Vector>
    class SemismoothNewton<Matrix, Vector, INVALID_BACKEND> final : public QPSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using Solver = utopia::LinearSolver<Matrix, Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;

    public:
        class ActiveTol {
        public:
            virtual ~ActiveTol() = default;

            ActiveTol(const Scalar tol) : tol_(tol) {}

            inline Scalar tol() const { return tol_; }

            Scalar &get() { return tol_; }

            inline void set(const Scalar tol) { tol_ = tol; }

            virtual Scalar tol(const SizeType i) const {
                UTOPIA_UNUSED(i);
                return tol_;
            }

            virtual void set_range(const Range &r) { UTOPIA_UNUSED(r); }

            virtual void became_active(const SizeType i) { UTOPIA_UNUSED(i); }

            virtual void became_inactive(const SizeType i) { UTOPIA_UNUSED(i); }

        private:
            Scalar tol_;
        };

        class AdaptiveActiveTol final : public ActiveTol {
        public:
            using ActiveTol::tol;

            AdaptiveActiveTol(const Scalar tol, const Range &r)
                : ActiveTol(tol), adaptive_tol_(r.extent(), tol), is_active_(r.extent(), false), r_begin_(r.begin()) {}

            AdaptiveActiveTol(const Scalar tol) : ActiveTol(tol), adaptive_tol_(1, tol), r_begin_(0) {}

            Scalar tol(const SizeType i) const override {
                assert(i < SizeType(adaptive_tol_.size()));
                return adaptive_tol_[i - r_begin_];
            }

            void set_range(const Range &r) override {
                r_begin_ = r.begin();
                adaptive_tol_.resize(r.extent());
                is_active_.resize(r.extent());
                clear();
            }

            // make it stickier everytime it becomes active
            void became_active(const SizeType i) override {
                auto idx = i - r_begin_;
                if (!is_active_[idx]) {
                    adaptive_tol_[idx] *= 2.0;
                    is_active_[idx] = true;
                }
            }

            void became_inactive(const SizeType i) override { is_active_[i - r_begin_] = false; }

            void clear() {
                std::fill(adaptive_tol_.begin(), adaptive_tol_.end(), tol());
                std::fill(is_active_.begin(), is_active_.end(), false);
            }

        private:
            std::vector<Scalar> adaptive_tol_;
            std::vector<bool> is_active_;
            SizeType r_begin_;
        };

        SemismoothNewton(const std::shared_ptr<Solver> &linear_solver)
            : linear_solver_(linear_solver),
              active_set_tol_(utopia::make_unique<ActiveTol>(1e-15)),
              linear_solve_zero_initial_guess_(true) {}

        SemismoothNewton *clone() const override {
            auto ptr = utopia::make_unique<SemismoothNewton>(std::shared_ptr<Solver>(linear_solver_->clone()));
            if (this->has_bound()) {
                ptr->set_box_constraints(this->get_box_constraints());
            }
            return ptr.release();
        }

        inline void set_active_set_tol(const Scalar tol) { active_set_tol_->set(tol); }

        void set_linear_solver(const std::shared_ptr<Solver> &ls) { linear_solver_ = ls; }

        void read(Input &in) override {
            QPSolver<Matrix, Vector>::read(in);
            bool use_adaptive_tol = false;
            in.get("use-adaptive-tol", use_adaptive_tol);

            if (use_adaptive_tol) {
                auto temp_tol = active_set_tol_->tol();
                active_set_tol_ = utopia::make_unique<AdaptiveActiveTol>(temp_tol);
            }

            in.get("active_set_tol", active_set_tol_->get());
            in.get("linear_solve_zero_initial_guess", linear_solve_zero_initial_guess_);
        }

        void print_usage(std::ostream &os) const override {
            QPSolver<Matrix, Vector>::print_usage(os);
            this->print_param_usage(os, "active_set_tol", "double", "Numerical tolerance.", "1e-15");
            this->print_param_usage(os,
                                    "linear_solve_zero_initial_guess",
                                    "bool",
                                    "Flag to reset initial guess for each linear solve.",
                                    "true");
        }

        inline bool smooth(const Vector &, Vector &) override { assert(false); return false; }

        inline bool apply(const Vector &b, Vector &x) override {
            const auto &A = *this->get_operator();

            assert(A.comm().size() == b.comm().size());
            assert(A.rows() == b.size());
            assert(x.empty() || (x.size() == b.size() && b.comm().size() == x.comm().size()));
            assert(this->has_bound());

            if (this->has_upper_bound() && this->has_lower_bound()) {
                assert(b.size() == this->upper_bound()->size());
                assert(b.size() == this->lower_bound()->size());

                box_solve(A, b, x);
            } else if ((this->has_upper_bound() && !this->has_lower_bound()) ||
                       (!this->has_upper_bound() && this->has_lower_bound())) {
                single_bound_solve(A, b, x);
            } else {
                std::cerr << "[Error] If you do not have constraints, use something else..... \n";
            }

            return true;
        }

        inline void init_memory(const Layout &layout) override {
            Super::init_memory(layout);
            if (linear_solver_) {
                linear_solver_->init_memory(layout);
            }
        }

        inline void update(const std::shared_ptr<const Matrix> &op) override {
            Super::update(op);

            if (linear_solver_) {
                linear_solver_->init_memory(row_layout(*op));
            }
        }

        inline void set_linear_solve_zero_initial_guess(const bool val) { linear_solve_zero_initial_guess_ = val; }

        inline bool linear_solve_zero_initial_guess() const { return linear_solve_zero_initial_guess_; }

    private:
        bool check_constraints(const Vector &x, const bool verbose) const {
            bool constrain_violated = false;
            if (this->has_upper_bound()) {
                const auto &upper_bound = this->get_upper_bound();
                Read<Vector> r_u(upper_bound);
                each_read(x, [&upper_bound, &constrain_violated, verbose](const SizeType i, const Scalar value) {
                    const Scalar upbo = upper_bound.get(i);
                    if (value > upbo + 1e-5) {
                        constrain_violated = true;
                        if (verbose) {
                            std::cerr << "[Warning] upper bound violated at: " << i << " " << value << " > " << upbo
                                      << std::endl;
                        }
                    }
                });
            }

            if (this->has_lower_bound()) {
                const auto &lower_bound = this->get_lower_bound();
                Read<Vector> r_u(lower_bound);
                each_read(x, [&lower_bound, &constrain_violated, verbose](const SizeType i, const Scalar value) {
                    const Scalar lobo = lower_bound.get(i);
                    if (value < lobo - 1e-5) {
                        constrain_violated = true;
                        if (verbose) {
                            std::cerr << "[Warning] lower bound violated at: " << i << " " << value << " < " << lobo
                                      << std::endl;
                        }
                    }
                });
            }

            return !constrain_violated;
        }

        static bool check_non_negative(const Vector &lambda, const bool verbose) {
            bool is_non_negative = true;
            each_read(lambda, [&is_non_negative, verbose](const SizeType i, const Scalar value) {
                if (value < -1e-5) {
                    is_non_negative = false;
                    if (verbose) {
                        std::cerr << "[Warning] negative value at: " << i << " " << value << std::endl;
                    }
                }
            });

            return is_non_negative;
        }

        static bool check_zero(const Vector &r) {
            const Scalar r_norm = norm2(r);

            if (r_norm > 1e-5) {
                std::cout << "non zero norm: " << r_norm << std::endl;
                return false;
            }

            return true;
        }

        static bool is_sane(const Matrix &A) {
            Vector d = diag(A);

            bool ret = true;

            each_read(d, [&ret](const SizeType /*i*/, const Scalar val) {
                if (std::abs(val) < 1e-16) {
                    ret = false;
                    std::cerr << "[Error] zero element on diagonal" << std::endl;
                }
            });

            return ret;
        }

        // We separate cases with 1 and 2 constraints in order to avoid usless computations in single constraint case
        bool single_bound_solve(const Matrix &A, const Vector &b, Vector &x_new) {
            bool is_upper_bound = this->has_upper_bound();

            assert(x_new.comm().size() == b.comm().size());
            assert(A.comm().size() == b.comm().size());

            Vector g;
            if (is_upper_bound) {
                g = this->get_upper_bound();

                assert(this->get_upper_bound().comm().size() == b.comm().size());
                assert(g.comm().size() == b.comm().size());
            } else {
                g = this->get_lower_bound();
            }

            active_set_tol_->set_range(range(g));

            auto A_layout = layout(A);
            auto v_layout = layout(b);

            // const SizeType local_N = local_size(x_new).get(0);

            SizeType iterations = 0;
            bool converged = false;

            Vector lambda(v_layout, 0.0);
            Vector active(v_layout, 0.0);

            Vector x_old = x_new;
            Vector d, prev_active;

            // active/inactive constraints
            Matrix A_c;
            Matrix I_c;

            if (is_sparse<Matrix>::value) {
                A_c.sparse(A_layout, 1, 1);
                I_c.sparse(A_layout, 1, 1);
            } else {
                A_c.dense(A_layout, 0.0);
                I_c.dense(A_layout, 0.0);
            }

            Scalar f_norm = 9e9;

            if (this->verbose()) this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "|| u_new - u_old ||"});

            while (!converged) {
                // reminder: complementarity-condition
                d = lambda + (x_new - g);
                {
                    Write<Vector> w_a(active);
                    Write<Matrix> w_A_c(A_c);
                    Write<Matrix> w_I_c(I_c);

                    Read<Vector> r_d(d);

                    const Range rr = row_range(A_c);

                    if (is_upper_bound) {
                        for (SizeType i = rr.begin(); i != rr.end(); i++) {
                            if (d.get(i) >= -active_set_tol_->tol(i)) {
                                A_c.set(i, i, 1.0);
                                active.set(i, 1.0);

                                active_set_tol_->became_active(i);

                                I_c.set(i, i, 0.0);
                            } else {
                                I_c.set(i, i, 1.0);
                                active.set(i, 0.0);

                                active_set_tol_->became_inactive(i);

                                A_c.set(i, i, 0.0);
                            }
                        }
                    } else {
                        // is_lower_bound
                        for (SizeType i = rr.begin(); i != rr.end(); i++) {
                            if (d.get(i) <= active_set_tol_->tol(i)) {
                                A_c.set(i, i, 1.0);
                                active.set(i, 1.0);

                                active_set_tol_->became_active(i);

                                I_c.set(i, i, 0.0);
                            } else {
                                I_c.set(i, i, 1.0);
                                active.set(i, 0.0);

                                active_set_tol_->became_inactive(i);

                                A_c.set(i, i, 0.0);
                            }
                        }
                    }
                }

                if (iterations > 0) {
                    const SizeType n_changed = Scalar(norm1(prev_active - active));

                    if (n_changed == 0) {
                        // active set doesn't change anymore => converged
                        // fix this to be done in other way
                        converged = this->check_convergence(iterations, 1e-15, 1, 1);
                        return true;
                    }
                }

                prev_active = active;

                Matrix H = A_c + I_c * A;
                Vector sub_g = (I_c * b + A_c * g);

                assert(g.comm().size() == x_new.comm().size());
                assert(sub_g.comm().size() == x_new.comm().size());

                assert(!has_nan_or_inf(H));
                assert(!has_nan_or_inf(g));
                assert(is_sane(H));

                if (this->linear_solve_zero_initial_guess()) {
                    x_new *= 0.;
                }

                if (!linear_solver_->solve(H, sub_g, x_new)) {
                    std::cerr << "[Error] linear solver did not manage to solve the linear system" << std::endl;
                    break;
                }

                if (has_nan_or_inf(x_new)) {
                    write("H.m", H);
                    write("s.m", sub_g);
                    write("x.m", x_new);

                    assert(!has_nan_or_inf(x_new));
                    std::cerr << "[Error] nan/inf entries in the solution vector" << std::endl;
                    converged = false;
                    break;
                }

                if (!check_zero(H * x_new - sub_g)) {
                    assert(check_zero(H * x_new - sub_g));
                    std::cerr << "[Error] non-zero residual returned by linear-solver" << std::endl;
                    converged = false;

                    write("H.m", H);
                    write("s.m", sub_g);
                    write("x.m", x_new);

                    break;
                }

                lambda = (b - A * x_new);
                assert(!has_nan_or_inf(lambda));

                f_norm = norm2(x_new - x_old);

                // print iteration status on every iteration
                if (this->verbose()) {
                    PrintInfo::print_iter_status(iterations, {f_norm});
                }

                converged = this->check_convergence(iterations, 1, 1, f_norm);

                x_old = x_new;
                iterations++;
            }

            // disp(lambda);
            assert(check_constraints(x_new, true));
            assert(check_non_negative(lambda, true));
            return converged;
        }

        bool box_solve(const Matrix &A, const Vector &b, Vector &x) {
            using namespace utopia;

            auto A_layout = layout(A);
            auto v_layout = row_layout(A);

            // const Size s_A = local_size(A);
            // const SizeType n = s_A.get(0);
            // const SizeType m = s_A.get(1);
            Scalar x_diff_norm = 0.0;

            SizeType it = 0;
            bool converged = false;

            const Vector &lb = this->get_lower_bound();
            const Vector &ub = this->get_upper_bound();

            active_set_tol_->set_range(range(ub));

            Vector lambda_p(v_layout, 0);
            Vector lambda_m(v_layout, 0);

            Vector active_m(v_layout, 0);
            Vector active_p(v_layout, 0);

            Vector x_old = x;
            Vector d_p, d_m;
            Vector prev_active_p, prev_active_m;
            Vector active;

            // active/inactive constraints
            Matrix A_c_p, A_c_m, A_s, I_c;
            Matrix H;

            if (this->verbose()) {
                this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "      || x_k - x_{k-1} ||"});
            }

            if (is_sparse<Matrix>::value) {
                A_c_p.identity(A_layout, 0.0);
                A_c_m.identity(A_layout, 0.0);
                I_c.identity(A_layout, 1.0);
            } else {
                A_c_p.dense(A_layout, 0.0);
                A_c_m.dense(A_layout, 0.0);
                I_c.identity(A_layout, 1.0);
            }

            while (!converged) {
                d_p = lambda_p + x - ub;
                d_m = lambda_m + x - lb;

                {
                    Write<Matrix> w_ac_p(A_c_p);
                    Write<Matrix> w_ac_m(A_c_m);
                    Write<Vector> w_a_m(active_m);
                    Write<Vector> w_a_p(active_p);

                    Read<Vector> rdp(d_p);
                    Read<Vector> rdm(d_m);
                    Read<Vector> rb(b);

                    Range rr = row_range(A_c_p);
                    for (SizeType i = rr.begin(); i != rr.end(); i++) {
                        if (d_p.get(i) >= -active_set_tol_->tol(i)) {
                            A_c_p.set(i, i, 1.0);

                            active_p.set(i, 1.0);
                            active_m.set(i, 0.0);

                            active_set_tol_->became_active(i);

                        } else if (d_m.get(i) <= active_set_tol_->tol(i)) {
                            A_c_m.set(i, i, 1.0);

                            active_m.set(i, 1.0);
                            active_p.set(i, 0.0);

                            active_set_tol_->became_active(i);

                        } else {
                            A_c_m.set(i, i, 0.0);
                            A_c_p.set(i, i, 0.0);

                            active_m.set(i, 0.0);
                            active_p.set(i, 0.0);

                            active_set_tol_->became_inactive(i);
                        }
                    }
                }

                active = active_m + active_p;
                A_s = diag(active);
                // I_c = local_identity(n, m);
                I_c *= 0.0;
                I_c.shift_diag(1.0);
                I_c -= A_s;

                if (it > 0) {
                    // active set doesn't change anymore => converged
                    const SizeType lb_changed = Scalar(norm1(prev_active_m - active_m));
                    const SizeType ub_changed = Scalar(norm1(prev_active_p - active_p));
                    const SizeType n_changed = lb_changed + ub_changed;

                    if (n_changed == 0) {
                        if (this->verbose() && mpi_world_rank() == 0) {
                            std::cout << "SemismoothNewton:: set is not changing ... \n";
                        }

                        converged = true;
                        break;
                    }
                }

                prev_active_m = active_m;
                prev_active_p = active_p;

                H = A_s + I_c * A;

                Vector rhs = I_c * b + A_c_p * ub + A_c_m * lb;
                linear_solver_->solve(H, rhs, x);

                lambda_p = A_c_p * (b - A * x);
                lambda_m = A_c_m * (b - A * x);

                // check for change in iterates
                x_diff_norm = norm2(x - x_old);

                // print iteration status on every iteration
                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {x_diff_norm});
                }

                if (!converged) {
                    converged = this->check_convergence(it, x_diff_norm, 1, 1);
                }

                x_old = x;
                it++;
            }

            assert(check_constraints(x, true));
            // assert(check_non_negative(lambda_p, true));
            return true;
        }

        std::shared_ptr<Solver> linear_solver_;
        std::unique_ptr<ActiveTol> active_set_tol_;
        bool linear_solve_zero_initial_guess_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_SEMISMOOTH_NEWTON_OLD_HPP
