#ifndef UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_IMPL_HPP
#define UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_IMPL_HPP

#include "utopia_PrimalInteriorPointSolver.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend>
    class PrimalInteriorPointSolver<Matrix, Vector, Backend>::Impl {
    public:
        void init(const Layout &layout) {
            correction.zeros(layout);
            buff.zeros(layout);

            if (linear_solver) {
                linear_solver->init_memory(layout);
            }
        }

        void pre_solve(const Vector &b, Vector &x) {
            if (box.has_upper_bound()) {
                x = min(x, *box.upper_bound());
            }

            // Lower bound is not supported yet though
            if (box.has_lower_bound()) {
                x = max(x, *box.lower_bound());
            }

            if (has_constraints_matrix()) {
                lambda.values(layout(x), 1.);
            } else {
                lambda = system_matrix() * x;
                lambda = b - lambda;
                lambda.e_max(min_val);
            }

            if (has_constraints_matrix()) {
                slack = *box.upper_bound() - constraints_matrix() * x;
            } else {
                slack = *box.upper_bound() - x;
            }

            slack.e_max(min_val);

            if (debug) {
                rename("s", slack);
                rename("l", lambda);
                rename("x", x);

                write("IP_pre_x.m", x);
                write("IP_pre_s.m", slack);
                write("IP_pre_l.m", lambda);
            }
        }

        void post_solve(Vector &x) {
            if (debug) {
                rename("s", slack);
                rename("l", lambda);
                rename("x", x);

                write("IP_post_x.m", x);
                write("IP_post_s.m", slack);
                write("IP_post_l.m", lambda);
            }
        }

        bool step(const Vector &b, Vector &x) {
            mu = duality_measure(lambda, slack);
            sigma = 0;

            residual_d(x, b, lambda, rd);
            residual_p(x, slack, rp);
            residual_s(lambda, slack, rs);

            if (empty(H)) {
                H = system_matrix();
            } else {
                H.same_nnz_pattern_copy(system_matrix());
            }

            e_pseudo_inv(slack, slack_inv);
            e_pseudo_inv(lambda, lambda_inv);

            Vector temp_vec = e_mul(slack_inv, lambda);

            ////////////////////////////////
            apply_selector(temp_vec);
            ////////////////////////////////

            if (has_constraints_matrix()) {
                // FIXME
                Matrix temp_mat = diag(temp_vec);
                Matrix temp_mat_2 = transpose(constraints_matrix()) * temp_mat * constraints_matrix();
                H += temp_mat_2;
            } else {
                H += diag(temp_vec);
            }

            if (has_constraints_matrix()) {
                // g = -transpose(constraints_matrix()) * e_mul(e_mul(rp, lambda) - rs, slack_inv);

                buff = e_mul(rp, lambda) - rs;
                buff = e_mul(buff, slack_inv);
                g = transpose(constraints_matrix()) * buff;
                g = -g;

            } else {
                g = e_mul(rp, lambda);
                g -= rs;
                g = e_mul(g, slack_inv);
                g = -g;
            }

            ////////////////////////////////
            apply_selector(g);
            ////////////////////////////////

            g -= rd;

            correction.set(0.);

            linear_solver->update(make_ref(H));

            if (!linear_solver->apply(g, correction)) {
                // break;
            }

            //////////////////////////////////////////////////////////////////////

            if (has_constraints_matrix()) {
                delta_lambda =
                    e_mul(slack_inv, e_mul(lambda, constraints_matrix() * correction) + e_mul(lambda, rp) - rs);
            } else {
                ////////////////////////////////
                apply_selector(correction);
                ////////////////////////////////

                delta_lambda = e_mul(lambda, correction);
                delta_lambda += e_mul(lambda, rp);
                delta_lambda -= rs;
                delta_lambda = e_mul(delta_lambda, slack_inv);
            }

            delta_slack = e_mul(slack, delta_lambda);
            delta_slack += rs;
            delta_slack = e_mul(delta_slack, lambda_inv);
            delta_slack = -delta_slack;

            // /////////////////////////////////////////////////////////////////////////
            // // https://en.wikipedia.org/wiki/Mehrotra_predictor%E2%80%93corrector_method

            Scalar alpha_primal = compute_alpha(slack, delta_slack);
            assert(alpha_primal == alpha_primal);

            Scalar alpha_dual = compute_alpha(lambda, delta_lambda);
            assert(alpha_dual == alpha_dual);

            Scalar alpha = std::min(alpha_primal, alpha_dual);
            buff = slack + alpha * delta_slack;
            buff = e_mul(buff, lambda + alpha * delta_lambda);
            Scalar mu_aff = sum(buff);

            mu_aff /= size(lambda);
            sigma = pow(mu_aff / mu, 3);
            assert(sigma == sigma);

            // //////////////////////////////////////
            rs_aff = e_mul(slack, lambda);
            buff = e_mul(delta_slack, delta_lambda);
            rs_aff += buff;
            rs_aff.shift(-sigma * mu);

            if (has_constraints_matrix()) {
                g = -transpose(constraints_matrix()) * e_mul(slack_inv, e_mul(lambda, rp) - rs_aff);
            } else {
                g = e_mul(lambda, rp) - rs_aff;
                g = e_mul(g, slack_inv);
                g = -g;
            }

            ////////////////////////////////
            apply_selector(g);
            ////////////////////////////////

            g -= rd;

            correction.set(0.);
            if (!linear_solver->apply(g, correction)) {
                // TODO check reason!
            }

            if (has_constraints_matrix()) {
                delta_lambda =
                    e_mul(slack_inv, e_mul(lambda, constraints_matrix() * correction) + e_mul(lambda, rp) - rs_aff);
            } else {
                delta_lambda = e_mul(lambda, correction);

                ////////////////////////////////
                apply_selector(delta_lambda);
                ////////////////////////////////

                buff = e_mul(lambda, rp);
                delta_lambda += buff;
                delta_lambda -= rs_aff;
                delta_lambda = e_mul(delta_lambda, slack_inv);
            }

            delta_slack = e_mul(slack, delta_lambda);
            delta_slack += rs_aff;
            delta_slack = e_mul(delta_slack, lambda_inv);
            delta_slack = -delta_slack;

            // ////////////////////////////////////////

            alpha_primal = compute_alpha_with_dumping(slack, delta_slack, dumping_parameter);
            assert(alpha_primal == alpha_primal);

            alpha_dual = compute_alpha_with_dumping(lambda, delta_lambda, dumping_parameter);
            assert(alpha_dual == alpha_dual);
            alpha = std::min(alpha_primal, alpha_dual);

            x += alpha * correction;
            lambda += alpha * delta_lambda;
            slack += alpha * delta_slack;

            // if (verbose) {
            //     std::stringstream ss;
            //     ss << "==============================\n";
            //     ss << "iteration:\t" << k << "\n";
            //     ss << "alpha:\t" << alpha << "\n";
            //     ss << "alpha_primal:\t" << alpha_primal << "\n";
            //     ss << "alpha_dual:\t" << alpha_dual << "\n";
            //     ss << "mu:\t" << mu << "\n";
            //     ss << "sigma:\t" << sigma << "\n";
            //     x.comm().root_print(ss.str());
            // }

            return true;
        }

        Scalar duality_measure(const Vector &lambda, const Vector &s) const { return dot(lambda, s) / size(lambda); }

        bool residual_d(const Vector &x, const Vector &b, const Vector &lambda, Vector &r_d) const {
            if (has_constraints_matrix()) {
                r_d = system_matrix() * x + transpose(constraints_matrix()) * lambda - b;
            } else if (has_selector()) {
                r_d = lambda;
                apply_selector(r_d);
                r_d += system_matrix() * x;
                r_d -= b;
            } else {
                r_d = system_matrix() * x;
                r_d += lambda;
                r_d -= b;
            }

            return true;
        }

        bool residual_p(const Vector &x, const Vector &s, Vector &r_p) const {
            if (has_constraints_matrix()) {
                r_p = constraints_matrix() * x - *box.upper_bound() + s;
            } else if (has_selector()) {
                r_p = x;
                apply_selector(r_p);
                r_p -= *box.upper_bound();
                r_p += s;
            } else {
                r_p = x - *box.upper_bound();
                r_p += s;
            }
            return true;
        }

        bool residual_s(const Vector &lambda, const Vector &s, Vector &r_s) {
            r_s = e_mul(s, lambda);

            if (sigma != 0) {
                r_s.shift(-sigma * mu);
            }

            return true;
        }

        bool is_positive(const Vector &x) const {
            auto x_view = local_view_device(x);

            int n_violations = 0;
            parallel_reduce(
                local_range_device(x),
                UTOPIA_LAMBDA(const SizeType i)->int {
                    auto xi = x_view.get(i);
                    return (xi <= 0) ? 1 : 0;
                },
                n_violations);

            n_violations = x.comm().sum(n_violations);
            return n_violations == 0;
        }

        Scalar compute_alpha_with_dumping(const Vector &x, const Vector &correction, const Scalar tau) {
            {
                auto x_view = local_view_device(x);
                auto correction_view = local_view_device(correction);

                buff.set(1.);

                auto buff_view = local_view_device(buff);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        auto xi = x_view.get(i);
                        auto dxi = correction_view.get(i);

                        if (xi + dxi < (1 - tau) * xi) {
                            auto val = -xi * tau / dxi;
                            buff_view.set(i, val);
                        }
                    });
            }

            return min(buff);
        }

        Scalar compute_alpha(const Vector &x, const Vector &correction) {
            {
                auto x_view = local_view_device(x);
                auto correction_view = local_view_device(correction);

                buff.set(1.);

                auto buff_view = local_view_device(buff);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        auto xi = x_view.get(i);
                        auto dxi = correction_view.get(i);

                        if (dxi < 0 && xi != 0) {
                            auto val = -xi / dxi;
                            buff_view.set(i, val);
                        }
                    });
            }

            return min(buff);
        }

        inline void apply_selector(Vector &vec) const {
            if (has_selector()) {
                vec = e_mul(*boolean_selector, vec);
            }
        }

        inline bool has_selector() const { return static_cast<bool>(boolean_selector); }

        inline const Matrix &system_matrix() const { return *system_matrix_; }

        inline bool has_constraints_matrix() const { return static_cast<bool>(constraints_matrix_); }
        inline const Matrix &constraints_matrix() const {
            assert(has_constraints_matrix());
            return *constraints_matrix_;
        }

        void determine_boolean_selector() const {
            this->box.determine_boolean_selector(-infinity, infinity, *boolean_selector);
        }

        std::shared_ptr<LinearSolver> linear_solver;

        Matrix H;
        Vector correction, residual;
        Vector lambda;
        Vector slack;
        Vector buff;

        Vector g;
        Vector delta_lambda, delta_slack;
        Vector rd, rp, rs, rs_aff;
        Vector slack_inv, lambda_inv;

        Scalar sigma{0};
        Scalar mu{0};
        Scalar dumping_parameter{1};
        Scalar min_val{1e-10};

        BoxConstraints<Vector> box;
        std::shared_ptr<const Matrix> constraints_matrix_;
        std::shared_ptr<const Matrix> system_matrix_;
        bool verbose{false};
        bool debug{false};

        std::shared_ptr<Vector> boolean_selector;
        bool auto_selector{false};
        Scalar infinity{10000};
    };

    template <class Matrix, class Vector, int Backend>
    PrimalInteriorPointSolver<Matrix, Vector, Backend>::PrimalInteriorPointSolver(
        std::shared_ptr<LinearSolver> linear_solver)
        : impl_(utopia::make_unique<Impl>()) {
        set_linear_solver(linear_solver);
    }

    template <class Matrix, class Vector, int Backend>
    PrimalInteriorPointSolver<Matrix, Vector, Backend>::PrimalInteriorPointSolver()
        : impl_(utopia::make_unique<Impl>()) {
        set_linear_solver(std::make_shared<OmniLinearSolver<Matrix, Vector>>());
    }

    template <class Matrix, class Vector, int Backend>
    PrimalInteriorPointSolver<Matrix, Vector, Backend>::~PrimalInteriorPointSolver() = default;

    template <class Matrix, class Vector, int Backend>
    PrimalInteriorPointSolver<Matrix, Vector, Backend> *PrimalInteriorPointSolver<Matrix, Vector, Backend>::clone()
        const {
        auto ptr = utopia::make_unique<PrimalInteriorPointSolver>(
            std::shared_ptr<LinearSolver>(impl_->linear_solver->clone()));
        if (this->has_bound()) {
            ptr->set_box_constraints(this->get_box_constraints());
        }
        return ptr.release();
    }

    template <class Matrix, class Vector, int Backend>
    void PrimalInteriorPointSolver<Matrix, Vector, Backend>::read(Input &in) {
        QPSolver<Matrix, Vector>::read(in);

        in.get("dumping_parameter", impl_->dumping_parameter);
        in.get("min_val", impl_->min_val);
        in.get("debug", impl_->debug);
        in.get("auto_selector", impl_->auto_selector);

        if (impl_->linear_solver) {
            in.get("linear_solver", *impl_->linear_solver);
        }
    }

    template <class Matrix, class Vector, int Backend>
    void PrimalInteriorPointSolver<Matrix, Vector, Backend>::print_usage(std::ostream &os) const {
        Super::print_usage(os);
    }

    template <class Matrix, class Vector, int Backend>
    bool PrimalInteriorPointSolver<Matrix, Vector, Backend>::smooth(const Vector &, Vector &) {
        assert(false && "IMPLEMENT ME");
        return false;
    }

    template <class Matrix, class Vector, int Backend>
    bool PrimalInteriorPointSolver<Matrix, Vector, Backend>::apply(const Vector &b, Vector &x) {
        assert(this->get_operator()->comm().size() == b.comm().size());
        assert(this->get_operator()->rows() == b.size());
        assert(x.empty() || (x.size() == b.size() && b.comm().size() == x.comm().size()));
        assert(this->has_upper_bound());

        impl_->box = this->get_box_constraints();

        if (impl_->auto_selector) {
            impl_->determine_boolean_selector();
        }

        if (this->verbose()) {
            this->init_solver("PrimalInteriorPointSolver comm.size = " + std::to_string(b.comm().size()),
                              {" it. ", "      || x_k - x_{k-1} ||", "duality_measure"});
        }

        impl_->verbose = this->verbose();

        impl_->pre_solve(b, x);

        SizeType it = 0;
        bool converged = false;
        while (!converged) {
            if (!impl_->step(b, x)) {
                break;
            }

            Scalar norm_c = norm2(impl_->correction);
            ++it;

            if (this->verbose()) {
                PrintInfo::print_iter_status(it, {norm_c, impl_->mu});
            }

            converged = this->check_convergence(it, 1, 1, norm_c);
        }

        impl_->post_solve(x);
        return converged;
    }

    template <class Matrix, class Vector, int Backend>
    void PrimalInteriorPointSolver<Matrix, Vector, Backend>::set_linear_solver(
        const std::shared_ptr<LinearSolver> &linear_solver) {
        impl_->linear_solver = linear_solver;
    }

    template <class Matrix, class Vector, int Backend>
    void PrimalInteriorPointSolver<Matrix, Vector, Backend>::init_memory(const Layout &layout) {
        Super::init_memory(layout);
        impl_->init(layout);
    }

    template <class Matrix, class Vector, int Backend>
    void PrimalInteriorPointSolver<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op) {
        Super::update(op);
        init_memory(row_layout(*op));

        // copy ptr of the matrix
        impl_->system_matrix_ = op;
    }

    template <class Matrix, class Vector, int Backend>
    void PrimalInteriorPointSolver<Matrix, Vector, Backend>::set_selection(
        const std::shared_ptr<Vector> &boolean_selector) {
        impl_->boolean_selector = boolean_selector;
    }

}  // namespace utopia

#endif  // UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_IMPL_HPP