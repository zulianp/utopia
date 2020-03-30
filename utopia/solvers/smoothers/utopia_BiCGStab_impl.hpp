#ifndef UTOPIA_BICG_STAB_IMPL_HPP
#define UTOPIA_BICG_STAB_IMPL_HPP

#include "utopia_BiCGStab.hpp"

namespace utopia {

    template<typename Matrix, typename Vector, int Backend>
    BiCGStab<Matrix, Vector, Backend>::BiCGStab(): initialized_(false)
    {}

    template<typename Matrix, typename Vector, int Backend>
    BiCGStab<Matrix, Vector, Backend> * BiCGStab<Matrix, Vector, Backend>::clone() const
    {
        return new BiCGStab(*this);
    }

    template<typename Matrix, typename Vector, int Backend>
    void BiCGStab<Matrix, Vector, Backend>::update(const Operator<Vector> &A)
    {
        Layout layout_rhs = row_layout(A);

        if(!initialized_ || !layout_rhs.same(layout_)) {
            init_memory(layout_rhs);
        }
    }


    template<typename Matrix, typename Vector, int Backend>
    void BiCGStab<Matrix, Vector, Backend>::init_memory(const Layout &layout)
    {
        OperatorBasedLinearSolver<Matrix, Vector>::init_memory(layout);

        v_.zeros(layout);
        p_.zeros(layout);
        h_.zeros(layout);
        y_.zeros(layout);
        t_.zeros(layout);
        s_.zeros(layout);
        z_.zeros(layout);
        K_inv_t_.zeros(layout);
    }

    //https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method (Preconditioned BiCGSTAB)
    template<typename Matrix, typename Vector, int Backend>
    bool BiCGStab<Matrix, Vector, Backend>::solve_preconditioned(const Operator<Vector> &A, const Vector &b, Vector &x)
    {
        const auto b_layout = layout(b);
        init_memory(b_layout);

        if(empty(x) || !b_layout.same(layout(x))) {
            x.zeros(b_layout);
        } else {
            assert( b_layout.same(layout(x)) );
        }

        Scalar rho = 1., rho_old = 1., alpha = 1., omega = 1., beta = 1.;
        A.apply(x, r_);
        r_ = b - r_;
        r0_ = r_;

        auto precond = this->get_preconditioner();

        this->init_solver("Utopia BiCGStab (preconditioned)", {"it. ", "||r||" });
        bool converged = false;

        SizeType it = 0;
        SizeType check_norm_each = 1;

        Scalar r_norm = 1e10;
        while(!converged)
        {
            rho = dot(r0_, r_);
            beta = (rho/rho_old) * (alpha/omega);
            p_ = r_ + beta * (p_ - omega * v_);
            precond->apply(p_, y_);
            A.apply(y_, v_);

            Scalar r0_dot_v = dot(r0_, v_);

            if(r0_dot_v == 0.0) {
                assert(false);

                A.apply(x, r_);
                r_ = b - r_;

                r_norm = norm2(s_);

                if(this->verbose()) {
                    PrintInfo::print_iter_status(it, { r_norm });
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
                break;
            }

            alpha = rho/r0_dot_v;

            h_ = x + alpha * y_;
            //if h is accurate enough
            //if(h_) { return }

            s_ = r_ - alpha * v_;

            if((it % check_norm_each) == 0) {
                r_norm = norm2(s_);

                if(this->verbose()) {
                    PrintInfo::print_iter_status(it, { r_norm });
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
                if(converged) { break; }
            }

            z_.set(0.0);
            precond->apply(s_, z_);
            A.apply(z_, t_);

            //FIXME should use left preconditioner instead
            precond->apply(t_, K_inv_t_);
            omega = dot(z_, K_inv_t_)/dot(K_inv_t_, K_inv_t_);

            assert(Scalar(dot(K_inv_t_, K_inv_t_)) != 0.0);

            if(std::isnan(omega)) {
                assert(false);
                break;
            }

            x = h_ + omega * z_;
            //if x is accurate enough
            //if(x_) { return }

            r_ = s_ - omega * t_;

            rho_old = rho;

            if((it % check_norm_each) == 0) {
                r_norm = norm2(r_);
                converged = this->check_convergence(it, r_norm, 1, 1);

                if(converged && this->verbose()) {
                    PrintInfo::print_iter_status(it, { r_norm });
                }
            }

            it++;
        }

        return converged;
    }

    //https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method (BiCGSTAB)
    template<typename Matrix, typename Vector, int Backend>
    bool BiCGStab<Matrix, Vector, Backend>::solve_unpreconditioned(const Operator<Vector> &A, const Vector &b, Vector &x)
    {
        const auto b_layout = layout(b);
        init_memory(b_layout);

        if(empty(x) || !b_layout.same(layout(x))) {
            x.zeros(b_layout);
        } else {
            assert( b_layout.same(layout(x)) );
        }

        Scalar rho = 1., rho_old = 1., alpha = 1., omega = 1., beta = 1.;
        A.apply(x, r_);
        r_ = b - r_;
        r0_ = r_;

        auto precond = this->get_preconditioner();

        this->init_solver("Utopia BiCGStab", {"it. ", "||r||" });
        bool converged = false;

        SizeType it = 0;
        SizeType check_norm_each = 1;

        Scalar r_norm = 1e10;
        while(!converged)
        {
            rho = dot(r0_, r_);

            if(rho == 0.) {
                converged = true;
                break;
            }

            beta = (rho/rho_old) * (alpha/omega);
            p_ = r_ + beta * (p_ - omega * v_);
            A.apply(p_, v_);

            const Scalar r0_dot_v = dot(r0_, v_);

            if(r0_dot_v == 0.0) {
                assert(false);
            }

            alpha = rho/r0_dot_v;

            h_ = x + alpha * p_;
                    //if h is accurate enough
                    //if(h_) { return }

            s_ = r_ - alpha * v_;

            // if((it % check_norm_each) == 0) {
            //     r_norm = norm2(s_);

            //     if(this->verbose()) {
            //         PrintInfo::print_iter_status({ Scalar(it), r_norm });
            //     }

            //     converged = this->check_convergence(it, r_norm, 1, 1);
            //     // if(converged) { break; }
            // }

            A.apply(s_, t_);
            omega = dot(t_, s_)/dot(t_, t_);

            if(std::isnan(omega) || std::isinf(omega)) {
                r_norm = norm2(s_);

                if(this->verbose()) {
                    PrintInfo::print_iter_status({ Scalar(it), r_norm });
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
                break;
            }

            x  = h_ + omega * s_;
            r_ = s_ - omega * t_;

            rho_old = rho;

            if((it % check_norm_each) == 0) {
                r_norm = norm2(r_);

                if(this->verbose()) {
                    PrintInfo::print_iter_status({ Scalar(it), r_norm });
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
            }

            it++;
        }

        return converged;
    }
}

#endif //UTOPIA_BICG_STAB_IMPL_HPP