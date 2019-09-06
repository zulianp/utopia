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
                            NewtonBase<Matrix, Vector>(linear_solver), tau_max_(1e14),reset_mass_matrix_(false), tau_zero_user_(-1)
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
        
        Scalar tau_init() const 
        { 
            return tau_zero_user_; 
        }
        
        void tau_init(const Scalar & tau_init)  
        {   
            tau_zero_user_ = tau_init; 
        }


        void reset_mass_matrix(const bool reset_mass)
        {
            reset_mass_matrix_ = reset_mass; 
        }


        void read(Input &in) override
        {
            NewtonBase<Matrix, Vector>::read(in);
            in.get("tau_max", tau_max_);
        }

        void print_usage(std::ostream &os) const override
        {
            NewtonBase<Matrix, Vector>::print_usage(os); 
            this->print_param_usage(os, "tau_max", "real", "Upper bound for tau.", "1e14"); 
        }

        void set_mass_matrix(const Matrix & M)
        {
            I_ = M; 
        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

            if(this->verbose())
                this->init_solver("PseudoContinuation", {" it. ", "|| g ||", "tau", "|| Delta x || "});

            bool converged = false; 
            SizeType it = 0; 

            Scalar g_norm, g_old, s_norm=9e9; 

            Vector g = local_zeros(local_size(x)), s; 
            Matrix H, H_damped; 

            fun.gradient(x, g); 
            g_norm = norm_l2(g); 

            fun.hessian(x, H); 

            if(empty(I_) || reset_mass_matrix_==true)
            {
                I_ = local_identity(local_size(H).get(0), local_size(H).get(1)); 
            }

            // tau = 1.0/g_norm; 
            Scalar tau = (tau_zero_user_ > 0)? tau_zero_user_ : std::max(1., 1./g_norm);

            // follows paper Combining TR methods and Rosenbrock Methods for Gradient systems
            // tau = std::min(g_norm, 10.0); 
            // tau = 1.0; 

            // option suitable to compare with AF-similarity solver
            // tau = 1.0/g_norm; 
            

            if(this->verbose())
                PrintInfo::print_iter_status(it, {g_norm, tau, 0.0});

            while(!converged)
            {

                H_damped = H + 1./tau * I_; 

                s = 0 * x; 
                this->linear_solve(H_damped, -1.0 * g, s);

                x += s; 
                
                fun.gradient(x, g); 

                g_old = g_norm; 
                // norms2(g, s, g_norm, s_norm); 
                s_norm = norm2(s); 
                g_norm = norm_l2(g); 
                
                if(g_norm > 1e-11)
                    tau = std::min(tau * g_old/g_norm, tau_max_); 
                else
                    tau = tau_max_; 

                it++; 

                if(this->verbose()){
                    PrintInfo::print_iter_status(it, {g_norm, tau, s_norm});
                }

                converged = this->check_convergence(it, g_norm, 9e9, s_norm);

                if(!converged)
                    fun.hessian(x, H); 

            } // outer solve loop while(!converged)

            return true;
        }
    
    private:
        Scalar norm_l2_2(const Vector & s)
        {   
            if(empty(I_)) {
                return dot(s, s);
            } else {
                return dot(s, I_*s); 
            }
        }

        Scalar norm_l2(const Vector & s)
        {   
            return std::sqrt(norm_l2_2(s)); 
        }

    private:
        Scalar tau_max_; 
        Matrix I_;
        bool reset_mass_matrix_; 
        Scalar tau_zero_user_; 


    };

}
#endif //UTOPIA_PSEUDO_CONTINUATION_HPP
