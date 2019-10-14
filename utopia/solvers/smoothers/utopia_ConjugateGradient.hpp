#ifndef UTOPIA_CONJUGATE_GRAD_H
#define UTOPIA_CONJUGATE_GRAD_H

#include "utopia_IterativeSolver.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Allocations.hpp"

#include <memory>

namespace utopia
{


    //FIXME also use the PreconditionedSolver interface properly
    /**
     * @brief      Conjugate Gradient solver. Works with all utopia tensor types.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class ConjugateGradient : public IterativeSolver<Matrix, Vector>, public Smoother<Matrix, Vector>, public MatrixFreeLinearSolver<Vector>
    {
        typedef UTOPIA_SCALAR(Vector) 	 Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::Preconditioner<Vector> Preconditioner;

    public:

        using IterativeSolver<Matrix, Vector>::solve;

        ConjugateGradient()
        : reset_initial_guess_(false), initialized_(false), loc_size_(0)
        {

        }

        void reset_initial_guess(const bool val)
        {
            reset_initial_guess_ = val;
        }


        void read(Input &in) override
        {
            IterativeSolver<Matrix, Vector>::read(in);
            Smoother<Matrix, Vector>::read(in);

            in.get("reset_initial_guess", reset_initial_guess_);

            if(precond_) {
                in.get("precond", *precond_);
            }
        }

        void print_usage(std::ostream &os) const override
        {
            IterativeSolver<Matrix, Vector>::print_usage(os);
            Smoother<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "reset_initial_guess", "bool", "Flag, which decides if initial guess should be reseted.", "false");
            this->print_param_usage(os, "precond", "Preconditioner", "Input parameters for preconditioner", "-");
        }


        /**
         * @brief      Solution routine for CG.
         *
         * @param[in]  b     The right hand side.
         * @param      x     The initial guess/solution.
         *
         * @return true if the linear system has been solved up to required tollerance. False otherwise
         */
        bool apply(const Vector &b, Vector &x) override
        {
            auto A_ptr = utopia::op(this->get_operator());

            SizeType loc_size_rhs = local_size(b); 
            if(!initialized_ || !b.comm().conjunction(loc_size_ == loc_size_rhs)) {
                    init(loc_size_rhs);
            }        

            return solve(*A_ptr, b, x);
        }


        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override
        {
            SizeType loc_size_rhs   = local_size(b); 
            if(!initialized_ || !b.comm().conjunction(loc_size_ == loc_size_rhs)) {
                    init(loc_size_rhs);
            }            

            if(precond_) {
                return preconditioned_solve(A, b, x);
            } else {
                return unpreconditioned_solve(A, b, x);
            }
        }

        /**
         * @brief      Sets the preconditioner.
         *
         * @param[in]  precond  The preconditioner.
         */
        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
        {
            precond_ = precond;
        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        void update(const std::shared_ptr<const Matrix> &op) override
        {
            IterativeSolver<Matrix, Vector>::update(op);

            if(precond_) {
                auto ls_ptr = dynamic_cast<LinearSolver<Matrix, Vector> *>(precond_.get());
                if(ls_ptr) {
                    ls_ptr->update(op);
                }
            }
        }

        bool smooth(const Vector &rhs, Vector &x) override
        {
            SizeType temp = this->max_it();
            this->max_it(this->sweeps());
            auto A_ptr = utopia::op(this->get_operator());
            unpreconditioned_solve(*A_ptr, rhs, x);
            this->max_it(temp);
            return true;
        }

        ConjugateGradient * clone() const override
        {
            return new ConjugateGradient(*this);
        }

    private:
        bool unpreconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x)
        {
            SizeType it = 0;
            Scalar rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9;

            assert(!empty(b));

            if(empty(x) || size(x) != size(b)) {
                x = local_zeros(local_size(b));
                r = b;
            } else {
                assert(local_size(x) == local_size(b));
                if(reset_initial_guess_) {
                    x.set(0.);
                }
                // r = b - A * x;
                UTOPIA_NO_ALLOC_BEGIN("CG:region1");
                A.apply(x, r);
                r = b - r;
                UTOPIA_NO_ALLOC_END();
            }

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

        bool preconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x)
        {
            SizeType it = 0;
            Scalar beta = 0., alpha = 1., r_norm = 9e9;

            z.set(0.0); 
            z_new.set(0.0); 

            if(empty(x) || size(x) != size(b)) {
                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region0");
                x = local_zeros(local_size(b));
                r = b;
                UTOPIA_NO_ALLOC_END();
            } else {
                assert(local_size(x) == local_size(b));

                if(reset_initial_guess_) {
                    x.set(0.);
                }

                UTOPIA_NO_ALLOC_BEGIN("CG_pre:region1");
                A.apply(x, r);
                r = b - r;
                UTOPIA_NO_ALLOC_END();
            }

            UTOPIA_NO_ALLOC_BEGIN("CG_pre:region2");
            precond_->apply(r, z);
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
                precond_->apply(r_new, z_new);
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

        bool check_solution(const Operator<Vector> &A, const Vector &x, const Vector &b) const
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

        void init(const SizeType &ls)
        {
            auto zero_expr = local_zeros(ls);

            //resets all buffers in case the size has changed
            r = zero_expr;
            p = zero_expr;
            q = zero_expr;
            Ap = zero_expr;
            r_new = zero_expr;
            z = zero_expr;
            z_new = zero_expr;

            initialized_ = true;    
            loc_size_ = ls;                    
        }

        std::shared_ptr<Preconditioner> precond_;
        Vector r, p, q, Ap, r_new, z, z_new;
        bool reset_initial_guess_;
        bool initialized_; 
        SizeType loc_size_;         
    };
}


#endif //UTOPIA_CONJUGATE_GRAD_H

