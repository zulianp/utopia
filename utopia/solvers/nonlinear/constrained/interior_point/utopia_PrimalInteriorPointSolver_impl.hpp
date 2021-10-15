#ifndef UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_IMPL_HPP
#define UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_IMPL_HPP

#include "utopia_PrimalInteriorPointSolver.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"

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

        void pre_solve(Vector &x) {
            x = min(x, *box.upper_bound());
            lambda.values(layout(x), 1.);

            if (has_constraints_matrix()) {
                slack = *box.upper_bound() - constraints_matrix() * x;
            } else {
                slack = *box.upper_bound() - x;
            }

            slack.e_max(1e-8);  // 0.1

            assert(is_positive(slack));
            assert(is_positive(lambda));
        }

        void post_solve() {}

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
                // H = system_matrix();
            }

            e_pseudo_inv(slack, slack_inv);
            e_pseudo_inv(lambda, lambda_inv);

            Vector temp_vec = e_mul(slack_inv, lambda);

            if (has_constraints_matrix()) {
                Matrix temp_mat = diag(temp_vec);
                Matrix temp_mat_2 = transpose(constraints_matrix()) * temp_mat * constraints_matrix();
                H += temp_mat_2;
            } else {
                H += diag(temp_vec);
            }

            if (has_constraints_matrix()) {
                g = -transpose(constraints_matrix()) * e_mul(e_mul(rp, lambda) - rs, slack_inv);
            } else {
                // g = -e_mul(e_mul(rp, lambda) - rs, slack_inv);

                g = e_mul(rp, lambda);
                g -= rs;
                g = e_mul(g, slack_inv);
                g = -g;
            }

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
                // delta_lambda = e_mul(slack_inv, e_mul(lambda, correction) + e_mul(lambda, rp) - rs);
                delta_lambda = e_mul(lambda, correction);
                delta_lambda += e_mul(lambda, rp);
                delta_lambda -= rs;
                delta_lambda = e_mul(delta_lambda, slack_inv);
            }

            // delta_slack = -e_mul(lambda_inv, rs + e_mul(slack, delta_lambda));

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

            // Scalar mu_aff = dot(slack + alpha * delta_slack, lambda + alpha * delta_lambda);

            buff = slack + alpha * delta_slack;
            buff = e_mul(buff, lambda + alpha * delta_lambda);
            Scalar mu_aff = sum(buff);

            mu_aff /= size(lambda);
            sigma = pow(mu_aff / mu, 3);
            assert(sigma == sigma);

            // //////////////////////////////////////
            // rs_aff = e_mul(slack, lambda) + e_mul(delta_slack, delta_lambda);
            rs_aff = e_mul(slack, lambda);
            buff = e_mul(delta_slack, delta_lambda);
            rs_aff += buff;
            rs_aff.shift(-sigma * mu);

            if (has_constraints_matrix()) {
                g = -transpose(constraints_matrix()) * e_mul(slack_inv, e_mul(lambda, rp) - rs_aff);
            } else {
                // g = -e_mul(slack_inv, e_mul(lambda, rp) - rs_aff);
                g = e_mul(lambda, rp) - rs_aff;
                g = e_mul(g, slack_inv);
                g = -g;
            }

            g -= rd;

            correction.set(0.);
            if (!linear_solver->apply(g, correction)) {
                // break;
            }

            if (has_constraints_matrix()) {
                delta_lambda =
                    e_mul(slack_inv, e_mul(lambda, constraints_matrix() * correction) + e_mul(lambda, rp) - rs_aff);
            } else {
                // delta_lambda = e_mul(slack_inv, e_mul(lambda, correction) + e_mul(lambda, rp) - rs_aff);
                delta_lambda = e_mul(lambda, correction);
                buff = e_mul(lambda, rp);
                delta_lambda += buff;
                delta_lambda -= rs_aff;
                delta_lambda = e_mul(delta_lambda, slack_inv);
            }

            // delta_slack = -e_mul(lambda_inv, rs_aff + e_mul(slack, delta_lambda));

            delta_slack = e_mul(slack, delta_lambda);
            delta_slack += rs_aff;
            delta_slack = e_mul(delta_slack, lambda_inv);
            delta_slack = -delta_slack;

            // ////////////////////////////////////////
            Scalar tau = 1;
            alpha_primal = compute_alpha_with_tau(slack, delta_slack, tau);
            assert(alpha_primal == alpha_primal);

            alpha_dual = compute_alpha_with_tau(lambda, delta_lambda, tau);
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
            } else {
                // r_d = system_matrix() * x + lambda - b;
                r_d = system_matrix() * x;
                r_d += lambda;
                r_d -= b;
            }

            return true;
        }

        bool residual_p(const Vector &x, const Vector &s, Vector &r_p) const {
            if (has_constraints_matrix()) {
                r_p = constraints_matrix() * x - *box.upper_bound() + s;
            } else {
                // r_p = x - *box.upper_bound() + s;
                r_p = x - *box.upper_bound();
                r_p += s;
            }
            return true;
        }

        bool residual_s(const Vector &lambda, const Vector &s, Vector &r_s) {
            r_s = e_mul(s, lambda);
            r_s.shift(-sigma * mu);
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

            // if (n_violations) utopia::out() << "n_violations: " << n_violations << "/" << x.size() << "\n";
            return n_violations == 0;
        }

        Scalar compute_alpha_with_tau(const Vector &x, const Vector &correction, const Scalar tau) {
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

        inline const Matrix &system_matrix() const { return *system_matrix_; }

        inline bool has_constraints_matrix() const { return static_cast<bool>(constraints_matrix_); }
        inline const Matrix &constraints_matrix() const {
            assert(has_constraints_matrix());
            return *constraints_matrix_;
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

        BoxConstraints<Vector> box;
        std::shared_ptr<const Matrix> constraints_matrix_;
        std::shared_ptr<const Matrix> system_matrix_;
        bool verbose{false};
    };

    template <class Matrix, class Vector, int Backend>
    PrimalInteriorPointSolver<Matrix, Vector, Backend>::PrimalInteriorPointSolver(
        std::shared_ptr<LinearSolver> linear_solver)
        : impl_(utopia::make_unique<Impl>()) {
        set_linear_solver(linear_solver);
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
        // const auto &A = *this->get_operator();

        assert(this->get_operator()->comm().size() == b.comm().size());
        assert(this->get_operator()->rows() == b.size());
        assert(x.empty() || (x.size() == b.size() && b.comm().size() == x.comm().size()));
        assert(this->has_upper_bound());
        // assert(!this->has_lower_bound());

        impl_->box = this->get_box_constraints();

        if (this->verbose()) {
            this->init_solver("PrimalInteriorPointSolver comm.size = " + std::to_string(b.comm().size()),
                              {" it. ", "      || x_k - x_{k-1} ||"});
        }

        impl_->verbose = this->verbose();

        impl_->pre_solve(x);

        SizeType it = 0;
        bool converged = false;
        while (!converged) {
            if (!impl_->step(b, x)) {
                break;
            }

            Scalar norm_c = norm2(impl_->correction);
            ++it;

            converged = this->check_convergence(it, 1, 1, norm_c);

            if (this->verbose()) {
                PrintInfo::print_iter_status(it, {norm_c});
            }
        }

        impl_->post_solve();
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

}  // namespace utopia

#endif  // UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_IMPL_HPP