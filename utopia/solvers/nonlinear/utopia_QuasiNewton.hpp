#ifndef UTOPIA_SOLVER_QUASI_NEWTON_HPP
#define UTOPIA_SOLVER_QUASI_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_QuasiNewtonBase.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Vector>
    class QuasiNewton : public QuasiNewtonBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                               Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                            SizeType;

        typedef utopia::HessianApproximation<Vector>                HessianApproximation;
        typedef utopia::MatrixFreeLinearSolver<Vector>              LinSolver;

    public:

        QuasiNewton(const std::shared_ptr <HessianApproximation> &hessian_approx,
                    const std::shared_ptr <LinSolver> &linear_solver):
                    QuasiNewtonBase<Vector>(hessian_approx, linear_solver), 
                    initialized_(false), 
                    loc_size_(0)
        {

        }

        bool solve(FunctionBase<Vector> &fun, Vector &x) override
        {
            using namespace utopia;

            Scalar g_norm, g0_norm, r_norm=1, s_norm=1;
            SizeType it = 0;
            Scalar alpha = 1.0;

            bool converged = false;

            SizeType loc_size_rhs = local_size(x); 
            if(!initialized_ || !g.comm().conjunction(loc_size_ == loc_size_rhs)) 
            {
                init_vectors(loc_size_rhs);
            }

            fun.gradient(x, g);
            g0_norm = norm2(g);
            g_norm = g0_norm;

            this->initialize_approximation(x, g); 

            

            if(this->verbose_) {
                this->init_solver("QUASI NEWTON", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm});
            }
            it++;            

            while(!converged)
            {
                UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:1");
                g_minus = -1.0 * g; 
                this->linear_solve(g_minus, s);
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:2");
                alpha = this->get_alpha(fun, g, x, s);
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:2.1");
                s *= alpha;
                x += s;
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:3");
                y = g;
                UTOPIA_NO_ALLOC_END();
                fun.gradient(x, g);

                UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:3.1");
                // norms needed for convergence check
                norms2(g, s, g_norm, s_norm);
                r_norm = g_norm/g0_norm;

                // diff between fresh and old grad...
                y = g - y;
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:4");
                this->update(s, y, x, g);
                UTOPIA_NO_ALLOC_END();

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm, alpha});

                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, r_norm, s_norm);

                it++;
            }

            this->print_statistics(it);
            return true;
        }



    private: 
        void init_vectors(const SizeType &ls)
        {
            auto zero_expr = local_zeros(ls);

            s       = zero_expr;
            g       = zero_expr;
            y       = zero_expr;
            g_minus = zero_expr;
                           
            initialized_ = true;    
            loc_size_ = ls;                                        
        }


        Vector g, s, y, g_minus;
        bool initialized_; 
        SizeType loc_size_;         

    };

}
#endif //UTOPIA_SOLVER_QUASI_NEWTON_HPP
