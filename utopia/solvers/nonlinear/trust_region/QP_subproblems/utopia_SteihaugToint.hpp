#ifndef UTOPIA_TR_L2_SUBPROBLEM_STEIHAUG_TOINT_HPP
#define UTOPIA_TR_L2_SUBPROBLEM_STEIHAUG_TOINT_HPP

#include "utopia_Allocations.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_TRSubproblem.hpp"

namespace utopia {
    /**
     * @brief      Class for Steihaug Toint conjugate gradient.
     */
    template <class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class SteihaugToint final : public OperatorBasedTRSubproblem<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        using OperatorBasedTRSubproblem<Matrix, Vector>::solve;
        using OperatorBasedTRSubproblem<Matrix, Vector>::update;

        using Preconditioner = utopia::Preconditioner<Vector>;

        SteihaugToint() : OperatorBasedTRSubproblem<Matrix, Vector>() {}

        void read(Input &in) override { OperatorBasedTRSubproblem<Matrix, Vector>::read(in); }

        void print_usage(std::ostream &os) const override {
            OperatorBasedTRSubproblem<Matrix, Vector>::print_usage(os);
        }

        SteihaugToint *clone() const override { return new SteihaugToint(*this); }

        void update(const Operator<Vector> &A) override {
            auto A_layout = row_layout(A);

            if (!initialized_ || !A_layout.same(layout_)) {
                init_memory(A_layout);
            }

            if (this->precond_) {
                auto ls_ptr = dynamic_cast<LinearSolver<Matrix, Vector> *>(this->precond_.get());
                if (ls_ptr) {
                    auto A_ptr = dynamic_cast<const Matrix *>(&A);
                    if (A_ptr) {
                        auto A_new = dynamic_cast<const Matrix &>(A);
                        ls_ptr->update(std::make_shared<const Matrix>(A_new));
                    }
                }
            }
        }

        bool apply(const Vector &b, Vector &x) override {
            minus_rhs = -1.0 * b;
            if (this->precond_) {
                // auto A_ptr = utopia::op(this->get_operator());
                return preconditioned_solve(*this->get_operator(), minus_rhs, x);
            } else {
                // auto A_ptr = utopia::op(this->get_operator());
                return unpreconditioned_solve(*this->get_operator(), minus_rhs, x);
            }
        }

        void use_precond_direction(const bool &use_precond_direction) {
            use_precond_direction_ = use_precond_direction;
        }

        bool use_precond_direction() { return use_precond_direction_; }

        // void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) {
        //   precond_ = precond;
        // }

        bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override {
            minus_rhs = rhs;
            minus_rhs *= -1.0;

            auto rhs_layout = layout(rhs);
            if (!initialized_ || !rhs_layout.same(layout_)) {
                init_memory(rhs_layout);
            }

            if (this->precond_) {
                return preconditioned_solve(A, minus_rhs, sol);
            } else {
                return unpreconditioned_solve(A, minus_rhs, sol);
            }
        }

    private:
        bool unpreconditioned_solve(const Operator<Vector> &B, const Vector &g, Vector &corr) {
            r = -1.0 * g;
            v_k = r;
            Scalar alpha, g_norm, d_B_d, z, z1;
            SizeType it = 1;

            corr.zeros(layout(g));
            g_norm = norm2(g);
            z = g_norm * g_norm;

            if (this->verbose()) {
                this->init_solver(" ST-CG ", {"it. ", "||r||"});
            }

            bool converged = false;

            // cudaProfilerStart();

            while (!converged) {
                B.apply(v_k, B_p_k);
                d_B_d = dot(v_k, B_p_k);

                if (d_B_d <= 0.0) {
                    Vector s = corr;
                    this->quad_solver(s, v_k, this->current_radius(), corr);
                    return true;
                }

                alpha = z / d_B_d;
                p_k = corr + alpha * v_k;

                if (norm2(p_k) >= this->current_radius()) {
                    Vector s = corr;
                    this->quad_solver(s, v_k, this->current_radius(), corr);
                    return true;
                }

                corr = p_k;
                r -= alpha * B_p_k;

                z1 = dot(r, r);
                v_k = r + (z1 / z) * v_k;

                z = z1;

                g_norm = std::sqrt(z);

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {g_norm});
                }

                converged = this->check_convergence(it, g_norm, 1, 1);
                it++;
            }

            // cudaProfilerStop();
            return true;
        }

        bool preconditioned_solve(const Operator<Vector> &B, const Vector &g, Vector &s_k) {
            bool converged = false;
            SizeType it = 0;

            if (empty(s_k))
                s_k.zeros(layout(g));
            else
                s_k.set(0.0);

            r = g;

            Scalar g_norm = 9e9;
            // cudaProfilerStart();

            this->init_solver(" Precond-ST-CG ", {"it. ", "||g||", "||s||", "||p||", "sMp"});
            if (this->verbose()) {
                PrintInfo::print_iter_status(it, {g_norm});
            }
            it++;

            if (empty(v_k))
                v_k.zeros(layout(g));
            else
                v_k.set(0.0);

            this->precond_->apply(r, v_k);

            p_k = -1.0 * v_k;

            Scalar alpha, kappa, betta;
            Scalar g_v_prod_old, g_v_prod_new;

            Scalar s_norm = 0.0, s_norm_new = 0.0, sMp = 0.0;
            Scalar r2 = this->current_radius() * this->current_radius();
            Scalar p_norm;

            if (use_precond_direction_) {
                if (this->compute_norm(it) || this->verbose() == true) {
                    dots(r, v_k, p_norm, r, r, g_norm);
                    g_norm = std::sqrt(g_norm);
                } else {
                    p_norm = dot(r, v_k);
                }
            } else {
                if (this->compute_norm(it) || this->verbose() == true) {
                    dots(p_k, p_k, p_norm, r, r, g_norm);
                    g_norm = std::sqrt(g_norm);
                } else {
                    p_norm = dot(p_k, p_k);
                }
            }

            // if preconditioner yields nans or inf, or is precond. dir is indefinite -
            // return gradient step if(!std::isfinite(p_norm) || p_norm < 0.0)
            if (!std::isfinite(p_norm)) {
                Scalar alpha_termination;
                if (r2 >= g_norm) {
                    alpha_termination = 1.0;  // grad. step is inside of tr boundary, just take it
                } else {
                    alpha_termination = std::sqrt(r2 / g_norm);  // grad. step is outside of
                                                                 // tr boundary, project on
                                                                 // the boundary
                }

                if (std::isfinite(alpha_termination)) {
                    s_k -= alpha_termination * r;
                }

                if (this->verbose() && mpi_world_rank() == 0) {
                    utopia::out() << "termination due to p_norm being nans/infs... \n";
                }

                this->check_convergence(it, g_norm, 1, 1e-15);

                return true;
            }

            while (!converged) {
                // UTOPIA_NO_ALLOC_BEGIN("STCG::region1");
                B.apply(p_k, B_p_k);
                // kappa = dot(p_k,B_p_k);

                dots(p_k, B_p_k, kappa, r, v_k, g_v_prod_old);

                // UTOPIA_NO_ALLOC_END();

                // identify negative curvature
                if (kappa <= 0.0) {
                    Scalar term1 = sMp * sMp + (p_norm * (r2 - s_norm));
                    Scalar tau = (std::sqrt(term1) - sMp) / p_norm;

                    if (std::isfinite(tau)) {
                        s_k += tau * p_k;
                    }

                    if (this->verbose() && mpi_world_rank() == 0) {
                        utopia::out() << "termination due to indefinite direction... \n";
                    }

                    this->check_convergence(it, g_norm, 1, 1e-15);

                    return true;
                }

                UTOPIA_NO_ALLOC_BEGIN("STCG::region2");
                // g_v_prod_old = dot(r, v_k);
                alpha = g_v_prod_old / kappa;

                s_norm_new = s_norm + (2.0 * alpha * sMp) + (alpha * alpha * p_norm);
                UTOPIA_NO_ALLOC_END();

                // ||s_k||_M > \Delta => terminate
                // norm squared should be used
                if (s_norm_new >= r2) {
                    UTOPIA_NO_ALLOC_BEGIN("STCG::region3");

                    Scalar term1 = sMp * sMp + (p_norm * (r2 - s_norm));
                    Scalar tau = (std::sqrt(term1) - sMp) / p_norm;

                    if (std::isfinite(tau)) {
                        s_k += tau * p_k;
                    }

                    if (this->verbose() && mpi_world_rank() == 0) {
                        utopia::out() << "termination due to correction exceeding TR radius... \n";
                    }

                    this->check_convergence(it, g_norm, 1, 1e-15);

                    UTOPIA_NO_ALLOC_END();
                    return true;
                }

                if (std::isfinite(alpha)) {
                    UTOPIA_NO_ALLOC_BEGIN("STCG::region4");
                    s_k += alpha * p_k;
                    UTOPIA_NO_ALLOC_END();
                } else {
                    return false;
                }

                UTOPIA_NO_ALLOC_BEGIN("STCG::region5");
                r += alpha * B_p_k;
                UTOPIA_NO_ALLOC_END();

                // apply preconditioner
                if (empty(v_k))
                    v_k.zeros(layout(r));
                else
                    v_k.set(0.0);

                UTOPIA_NO_ALLOC_BEGIN("STCG::region6");
                this->precond_->apply(r, v_k);
                UTOPIA_NO_ALLOC_END();

                g_v_prod_new = dot(r, v_k);

                // if preconditioner yields nans or inf, or is precond. dir is indefinite
                // - just return current step
                if (!std::isfinite(g_v_prod_new)) {
                    return true;
                }

                betta = g_v_prod_new / g_v_prod_old;

                UTOPIA_NO_ALLOC_BEGIN("STCG::region6.2");
                p_k = betta * p_k - v_k;
                UTOPIA_NO_ALLOC_END();

                // updating norms recursively  - see TR book
                if (use_precond_direction_) {
                    sMp = betta * (sMp + (alpha * p_norm));
                    p_norm = g_v_prod_new + (betta * betta * p_norm);
                    s_norm = s_norm_new;
                } else {
                    UTOPIA_NO_ALLOC_BEGIN("STCG::region7");
                    dots(p_k, s_k, sMp, p_k, p_k, p_norm, s_k, s_k, s_norm);
                    UTOPIA_NO_ALLOC_END();
                }

                // TODO:: check if there is something else possible
                if (this->compute_norm(it) || this->verbose() == true) {
                    g_norm = norm2(r);
                }

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, p_norm, sMp});
                }

                if (!std::isfinite(g_norm)) {
                    return false;
                }

                converged = this->check_convergence(it, g_norm, 1, 1);
                it++;
            }

            // cudaProfilerStop();
            return true;
        }

    public:
        void init_memory(const Layout &layout) override {
            // resets all buffers in case the size has changed
            v_k.zeros(layout);
            r.zeros(layout);
            p_k.zeros(layout);
            B_p_k.zeros(layout);
            minus_rhs.zeros(layout);

            initialized_ = true;
            layout_ = layout;
        }

    private:
        Vector v_k, r, p_k, B_p_k, minus_rhs;
        // std::shared_ptr<Preconditioner> precond_; /*!< Preconditioner to be used. */
        bool use_precond_direction_{false};
        bool initialized_{false};
        Layout layout_;
    };

}  // namespace utopia

#endif  // UTOPIA_TR_SUBPROBLEM_STEIHAUG_TOINT_HPP
