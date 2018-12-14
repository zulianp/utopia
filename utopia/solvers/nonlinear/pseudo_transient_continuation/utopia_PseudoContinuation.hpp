#ifndef UTOPIA_PSEUDO_CONTINUATION_HPP
#define UTOPIA_PSEUDO_CONTINUATION_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{   
    /**
     * @brief This solver is implementation of following papers: 
     *  Convergence analysis of pseudo-transient contiuation by Kelley, Keyes
     *  Pseudotransient Continuation and Differential-Algebraic equations by Colley, Kelley, Keyes
     *  Projected pseudo-transient continuation by Kelley, Liao, Qi, Chu, Reese, Winton
     *  
     */   
    template<class Matrix, class Vector>
    class PseudoContinuation final: public NewtonBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef typename NewtonBase<Matrix, Vector>::Solver Solver;

        using NewtonBase<Matrix, Vector>::print_statistics; 


    public:
       PseudoContinuation(  const std::shared_ptr <Solver> &linear_solver):
                            NewtonBase<Matrix, Vector>(linear_solver), tau_max_(1e14)
        {

        }

        void tau_max(const Scalar & tau_max)
        {
            tau_max_ = tau_max; 
        }

        Scalar tau_max() const
        {
            return tau_max_; 
        }

        void read(Input &in) override
        {
            NewtonBase<Matrix, Vector>::read(in);
            in.get("tau_max", tau_max_);
        }


        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

            if(this->verbose())
                this->init_solver("PseudoContinuation", {" it. ", "|| g ||", "tau", "|| Delta x || "});

            bool converged = false; 
            SizeType it = 0; 

            Scalar g_norm, g_old, s_norm=9e9, tau; 

            Vector g = local_zeros(local_size(x)), s; 
            Matrix H, H_damped; 

            fun.gradient(x, g); 
            g_norm = norm2(g); 

            fun.hessian(x, H); 
            Matrix I = local_identity(local_size(H)); 

            // tau = 1.0/g_norm; 
            tau = g_norm; 

            if(this->verbose())
                PrintInfo::print_iter_status(it, {g_norm, tau, 0.0});

            while(!converged)
            {
                H_damped = H + 1./tau * I; 

                s = 0 * x; 
                this->linear_solve(H_damped, -1.0 * g, s);
                x += s; 

                fun.gradient(x, g); 

                g_old = g_norm; 
                norms2(g, s, g_norm, s_norm); 
                tau = std::min(tau * g_old/g_norm, tau_max_); 

                it++; 

                if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm, tau, s_norm});

                converged = this->check_convergence(it, g_norm, 9e9, s_norm);

                if(!converged)
                    fun.hessian(x, H); 

            } // outer solve loop while(!converged)

            return true;
        }
    
    private:
        Scalar tau_max_; 


    };

}
#endif //UTOPIA_PSEUDO_CONTINUATION_HPP
