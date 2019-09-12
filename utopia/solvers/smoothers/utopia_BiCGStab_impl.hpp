#ifndef UTOPIA_BICG_STAB_IMPL_HPP
#define UTOPIA_BICG_STAB_IMPL_HPP

#include "utopia_BiCGStab.hpp"

namespace utopia {

    template<typename Matrix, typename Vector, int Backend>
    BiCGStab<Matrix, Vector, Backend>::BiCGStab()
    {}

    template<typename Matrix, typename Vector, int Backend>
    BiCGStab<Matrix, Vector, Backend> * BiCGStab<Matrix, Vector, Backend>::clone() const
    {
        return new BiCGStab(*this);
    }

    template<typename Matrix, typename Vector, int Backend>
    void BiCGStab<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op)
    {
        PreconditionedSolver<Matrix, Vector>::update(op);
        // init({local_size(*op).get(0)});
    }

    template<typename Matrix, typename Vector, int Backend>
    void BiCGStab<Matrix, Vector, Backend>::init(const Size &ls)
    {
        v_ = local_zeros(ls);
        p_ = local_zeros(ls);
        h_ = local_zeros(ls);
        y_ = local_zeros(ls);
        t_ = local_zeros(ls);
        s_ = local_zeros(ls);
        z_ = local_zeros(ls);
        K_inv_t_ = local_zeros(ls);
    }

    template<typename Matrix, typename Vector, int Backend>
    bool BiCGStab<Matrix, Vector, Backend>::smooth(const Vector &rhs, Vector &x)
    {
        SizeType temp = this->max_it();
        this->max_it(this->sweeps());
        auto A_ptr = utopia::op(this->get_operator());
        solve_unpreconditioned(*A_ptr, rhs, x);
        this->max_it(temp);
        return true;
    }

    //https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method (Preconditioned BiCGSTAB)
    template<typename Matrix, typename Vector, int Backend>
    bool BiCGStab<Matrix, Vector, Backend>::solve_preconditioned(const Operator<Vector> &A, const Vector &b, Vector &x)
    {
        const auto ls = local_size(b);
        init(ls);

        if(empty(x) || size(x) != size(b)) {
            x = local_zeros(ls);
        } else {
            assert(local_size(x) == ls);
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
            }

            alpha = rho/r0_dot_v;

            h_ = x + alpha * y_;
            //if h is accurate enough
            //if(h_) { return }

            s_ = r_ - alpha * v_;

            if((it % check_norm_each) == 0) {
                r_norm = norm2(s_);

                if(this->verbose()) {
                    PrintInfo::print_iter_status({ Scalar(it), r_norm });
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
                if(converged) { break; }
            }

            precond->apply(s_, z_);
            A.apply(z_, t_);

            //FIXME should use left preconditioner instead
            precond->apply(t_, K_inv_t_);
            omega = dot(z_, K_inv_t_)/dot(K_inv_t_, K_inv_t_);

            assert(Scalar(dot(K_inv_t_, K_inv_t_)) != 0.0);

            if(std::isnan(omega)) {
                assert(false);
            }

            x = h_ + omega * z_;
            //if x is accurate enough
            //if(x_) { return }

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

    //https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method (BiCGSTAB)
    template<typename Matrix, typename Vector, int Backend>
    bool BiCGStab<Matrix, Vector, Backend>::solve_unpreconditioned(const Operator<Vector> &A, const Vector &b, Vector &x)
    {
        const auto ls = local_size(b);
        init(ls);

        if(empty(x) || size(x) != size(b)) {
            x = local_zeros(ls);
        } else {
            assert(local_size(x) == ls);
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

            if((it % check_norm_each) == 0) {
                r_norm = norm2(s_);

                if(this->verbose()) {
                    PrintInfo::print_iter_status({ Scalar(it), r_norm });
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
                if(converged) { break; }
            }

            A.apply(s_, t_);
            omega = dot(t_, s_)/dot(t_, t_);

            if(std::isnan(omega) || std::isinf(omega)) {
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