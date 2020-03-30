#ifndef UTOPIA_CONJUGATE_GRADIENT_IMPL_HPP
#define UTOPIA_CONJUGATE_GRADIENT_IMPL_HPP

#include "utopia_ConjugateGradient.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Allocations.hpp"

#include <memory>

namespace utopia {

    template<class Matrix, class Vector, int Backend>
    ConjugateGradient<Matrix, Vector, Backend>::ConjugateGradient()
    : reset_initial_guess_(false), initialized_(false)
    {}

    template<class Matrix, class Vector, int Backend>
    void ConjugateGradient<Matrix, Vector, Backend>::reset_initial_guess(const bool val)
    {
        reset_initial_guess_ = val;
    }

    template<class Matrix, class Vector, int Backend>
    void ConjugateGradient<Matrix, Vector, Backend>::read(Input &in)
    {
        OperatorBasedLinearSolver<Matrix, Vector>::read(in);
        in.get("reset_initial_guess", reset_initial_guess_);
    }

    template<class Matrix, class Vector, int Backend>
    void ConjugateGradient<Matrix, Vector, Backend>::print_usage(std::ostream &os) const
    {
        OperatorBasedLinearSolver<Matrix, Vector>::print_usage(os);
        this->print_param_usage(os, "reset_initial_guess", "bool", "Flag, which decides if initial guess should be reseted.", "false");
    }

    template<class Matrix, class Vector, int Backend>
    bool ConjugateGradient<Matrix, Vector, Backend>::solve(const Operator<Vector> &A, const Vector &b, Vector &x)
    {
        if(this->has_preconditioner()) {
            return preconditioned_solve(A, b, x);
        } else {
            return unpreconditioned_solve(A, b, x);
        }
    }

    template<class Matrix, class Vector, int Backend>
    void ConjugateGradient<Matrix, Vector, Backend>::update(const Operator<Vector> &A)
    {
        const auto layout_rhs = row_layout(A);

        if(!initialized_ || !layout_rhs.same(layout_)) {
            init_memory(layout_rhs);
        }
    }

    template<class Matrix, class Vector, int Backend>
    ConjugateGradient<Matrix, Vector, Backend> * ConjugateGradient<Matrix, Vector, Backend>::clone() const
    {
        return new ConjugateGradient(*this);
    }

    template<class Matrix, class Vector, int Backend>
    void ConjugateGradient<Matrix, Vector, Backend>::gradient_descent_step(
        const Operator<Vector> &A,
        const Vector &b,
        Vector &x)
    {
        A.apply(x, r);
        //-grad = b - A * x
        r = b - r;
        x += r;
    }

    template<class Matrix, class Vector, int Backend>
    bool ConjugateGradient<Matrix, Vector, Backend>::unpreconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x)
    {
        SizeType it = 0;
        Scalar rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9;

        assert(!empty(b));

        if(empty(x) || size(x) != size(b)) {
            x = local_zeros(local_size(b));
        } else {
            assert(local_size(x) == local_size(b));
            if(reset_initial_guess_) {
                x.set(0.);
            }
        }

        gradient_descent_step(A, b, x);

        // r = b - A * x;
        UTOPIA_NO_ALLOC_BEGIN("CG:region1");
        A.apply(x, r);
        r = b - r;
        UTOPIA_NO_ALLOC_END();
        // }

        this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||" });
        bool converged = false;

        SizeType check_norm_each = 1;

        while(!converged)
        {
            rho = dot(r, r);

            if(rho == 0.) {
                converged = true;
                break;
            }

            if(it > 0)
            {
                beta = rho/rho_1;
                UTOPIA_NO_ALLOC_BEGIN("CG:region2");
                p = r + beta * p;
                UTOPIA_NO_ALLOC_END();
            }
            else
            {
                UTOPIA_NO_ALLOC_BEGIN("CG:region3");
                p = r;
                UTOPIA_NO_ALLOC_END();
            }

            UTOPIA_NO_ALLOC_BEGIN("CG:region4");
            // q = A * p;
            A.apply(p, q);

            Scalar dot_pq = dot(p, q);

            UTOPIA_NO_ALLOC_END();

            if(dot_pq == 0.) {
                //TODO handle properly
                utopia_warning("prevented division by zero");
                converged = true;
                break;
            }

            UTOPIA_NO_ALLOC_BEGIN("CG:region5");
            alpha = rho / dot_pq;

            x += alpha * p;
            r -= alpha * q;

            rho_1 = rho;
            UTOPIA_NO_ALLOC_END();

            if((it % check_norm_each) == 0) {
                // r =
                // A.apply(x, r);
                // r = b - r;

                r_norm = norm2(r);

                if(this->verbose()) {
                    PrintInfo::print_iter_status(it, { r_norm });
                }

                converged = this->check_convergence(it, r_norm, 1, 1);
            }

            it++;
        }

        return converged;
    }

    template<class Matrix, class Vector, int Backend>
    bool ConjugateGradient<Matrix, Vector, Backend>::preconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x)
    {
        SizeType it = 0;
        Scalar beta = 0., alpha = 1., r_norm = 9e9;

        z.set(0.0);
        z_new.set(0.0);

        auto precond = this->get_preconditioner();

        if(empty(x) || size(x) != size(b)) {
            UTOPIA_NO_ALLOC_BEGIN("CG_pre:region0");
            x = local_zeros(local_size(b));
            // r = b;
            // UTOPIA_NO_ALLOC_END();
        } else {
            assert(local_size(x) == local_size(b));

            if(reset_initial_guess_) {
                x.set(0.);
            }
        }

        gradient_descent_step(A, b, x);

        UTOPIA_NO_ALLOC_BEGIN("CG_pre:region1");
        A.apply(x, r);
        r = b - r;
        UTOPIA_NO_ALLOC_END();


        UTOPIA_NO_ALLOC_BEGIN("CG_pre:region2");
        precond->apply(r, z);
        p = z;
        UTOPIA_NO_ALLOC_END();

        this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||" });
        bool stop = false;

        while(!stop)
        {
            // Ap = A*p;
            UTOPIA_NO_ALLOC_BEGIN("CG_pre:region3");
            A.apply(p, Ap);
            alpha = dot(r, z)/dot(p, Ap);
            UTOPIA_NO_ALLOC_END();

            if(std::isinf(alpha) || std::isnan(alpha)) {
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

            if(r_norm < this->atol()) {
                if(this->verbose()) {
                    PrintInfo::print_iter_status(it, {r_norm});
                }

                stop = this->check_convergence(it, r_norm, 1, 1);
                break;
            }

            UTOPIA_NO_ALLOC_BEGIN("CG_pre:region5");
            z_new.set(0.0);
            precond->apply(r_new, z_new);
            beta = dot(z_new, r_new)/dot(z, r);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("CG_pre:region5.1");
            p = z_new + beta * p;
            r = r_new;
            z = z_new;
            UTOPIA_NO_ALLOC_END();

            if(this->verbose()) {
                PrintInfo::print_iter_status(it, {r_norm});
            }

            stop = this->check_convergence(it, r_norm, 1, 1);
            it++;
        }

        if(r_norm <= this->atol()) {
            //FIXME sometimes this fails for some reason
            // assert(check_solution(A, x, b));
            return true;
        } else {
            return false;
        }
    }

    template<class Matrix, class Vector, int Backend>
    bool ConjugateGradient<Matrix, Vector, Backend>::check_solution(const Operator<Vector> &A, const Vector &x, const Vector &b) const
    {
        Vector r;
        A.apply(x, r);
        r -= b;

        const Scalar r_norm = norm2(r);

        if(r_norm > 100 * this->atol()) {
            // write("A.m", *this->get_operator());
            // disp(*this->get_operator());
            assert(r_norm <= this->atol());
            return false;
        }

        return true;
    }

    template<class Matrix, class Vector, int Backend>
    void ConjugateGradient<Matrix, Vector, Backend>::init_memory(const Layout &layout)
    {
        assert(layout.local_size() > 0);
        OperatorBasedLinearSolver<Matrix, Vector>::init_memory(layout);

        //resets all buffers in case the size has changed
        r.zeros(layout);
        p.zeros(layout);
        q.zeros(layout);
        Ap.zeros(layout);
        r_new.zeros(layout);
        z.zeros(layout);
        z_new.zeros(layout);

        initialized_ = true;
        layout_ = layout;
    }

}


#endif //UTOPIA_CONJUGATE_GRADIENT_IMPL_HPP
