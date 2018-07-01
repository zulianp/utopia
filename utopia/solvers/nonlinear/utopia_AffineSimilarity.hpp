#ifndef UTOPIA_AFFINE_SIMILARITY_TRANSFORM_HPP
#define UTOPIA_AFFINE_SIMILARITY_TRANSFORM_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class AffineSimilarity : public NonLinearSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver Solver;
        typedef utopia::LSStrategy<Matrix, Vector> LSStrategy; 

    public:
       AffineSimilarity(    const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(), 
                            const Parameters params                       = Parameters() ):
                            NonLinearSolver<Matrix, Vector>(linear_solver, params), mass_init_(false)
                            {
                                //set_parameters(params);
                            }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

            Vector g, s, rhs, x_trial;
            Matrix H, A;

            Scalar g_norm, g0_norm=9e9, s_norm=9e9, tau;
            SizeType it = 0;

            bool converged = false, taken=1;

            fun.gradient(x, g);
            g = -1.0*g;

            g0_norm = norm2(g);
            g_norm = g0_norm;

            this->init_solver("Affine similarity", {" it. ", "|| g ||", "|| s_k || ", "tau", "nu",  "taken"});
            tau = 1.0/norm2(g); 
            
            if(this->verbose_)
                PrintInfo::print_iter_status(it, {g_norm, 0, tau});
            it++;

            Vector g_old = g; 

            while(!converged)
            {
                // if(taken == 1)
                // {
                    fun.hessian(x, H);
                    H *= -1.0; 

                    fun.gradient(x, g);                
                    g *= -1.0;       
                    g0_norm = norm2(g);          
                // }
                
                A = H - 1.0/tau * M_; 
                rhs = -1.0 * g; 
                
                //find direction step
                s = local_zeros(local_size(x));
                this->linear_solve(A, rhs, s);

                x_trial = x+s; 

                // plain correction, without step-size
                s = 1.0/tau * s; 
                s_norm = norm2(s); 
                
                
                Scalar nu = dot(s, H* s)/(s_norm*s_norm*tau);

                // gradient of x_trial 
                Vector g_trial; 
                fun.gradient(x_trial, g_trial);  
                g_trial *= -1.0;     

                // difference between gradient of trial point and correction
                Vector gs_diff = (g_trial - (M_ * s)); 
                Scalar nom = dot(s, ( (1.0/tau * M_ * s) - g)); 

                Scalar help_denom = (2.0 * norm2(gs_diff) * s_norm); 
                tau = tau  *  std::abs(nom)/ help_denom; 
 
                
                {
                    if(norm2(g_trial) < norm2(g))
                    {
                        x = x_trial; 
                        taken = 1; 
                    }
                    else
                    {
                        taken=0;
                    }


                    // clamping values of tau
                    if(std::isinf(tau))
                        tau = 9e249;        
                    else if(std::isnan(tau))
                        tau = 1e-10;
                    else if(tau==0)
                        tau = 1e-10;                                        
                }
 
            
                fun.gradient(x, g);  
                g *= -1.0;     
                g_norm = norm2(g);

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, tau, nu,  double(taken)});

                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                it++;
            }

            return true;
        }


        void set_mass_matrix(const Matrix & M)
        {
            M_ = M; 
            mass_init_ = true; 
        }


    private:
        Matrix M_;  // mass matrix 
        bool mass_init_; 



    };

}
#endif //UTOPIA_AFFINE_SIMILARITY_TRANSFORM_HPP
