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

            Scalar g_norm, g0_norm=9e9, r_norm=9e9, s_norm=9e9, tau;
            SizeType it = 0;

            bool converged = false, taken=0;

            fun.gradient(x, g);
            g = -1.0*g;

            g0_norm = norm2(g);
            g_norm = g0_norm;


            this->init_solver("Affine similarity", {" it. ", "|| g ||", "<s,g>", "|| s_k || ", "tau", "nu",  "taken"});

            tau = norm2(g); 
            tau = 1/tau; 
            
            if(this->verbose_)
                PrintInfo::print_iter_status(it, {g_norm, 1, 0, tau});
            it++;

            SizeType inner_counter = 0; 

            fun.hessian(x, H);
            H = -1.0* H; 

            Scalar nu_zero_approx = dot(g, H*g)/(g0_norm*g0_norm); 
            std::cout<<"nu_zero_approx: "<< nu_zero_approx << "  \n"; 

            Vector g_old = g; 

            while(!converged)
            {
                if(taken == 1)
                {
                    fun.hessian(x, H);
                    H *= -1.0; 

                    fun.gradient(x, g);                
                    g *= -1.0;       
                    g0_norm = norm2(g);          
                }
                
                A = M_ - tau*H; 
                rhs = tau*g; 
                
                //find direction step
                s = local_zeros(local_size(x));
                this->linear_solve(A, rhs, s);

                x_trial = x+s; 

                // plain correction, without step-size
                Vector t0 = (H*s); 
                Scalar inv_tau = 1.0/tau; 
                s = inv_tau * s; 
                s_norm = norm2(s); 

                Scalar s_test_norm = norm2(M_ * s); 

                if(s_test_norm > norm2(g))
                    std::cout<<"iteration should terminate.... \n"; 

                
                Scalar nu = dot(s, ( (M_ * s) - g))/(s_norm*s_norm*tau);


                // gradient of x_trial 
                fun.gradient(x_trial, g);  
                g *= -1.0;     

                // difference between gradient of trial point and correction
                Vector gs_diff = g - (M_ *s); 
                Scalar L2 = 2.0 * norm2(gs_diff); 


                Scalar help = (tau*tau * s_norm*s_norm); 
                L2 = L2/ help;
                
                Scalar help2 = L2 * s_norm; 
                tau = std::abs(nu)/help2; 

                Scalar denom = (g0_norm * L2);
                denom = denom - (nu*nu); 
                std::cout<<"denom: "<< denom << "  \n"; 
                Scalar tau_opt = std::abs(nu)/ denom; 

                std::cout<<"tau_opt: "<< tau_opt << "  tau: "<< tau << "   \n"; 

                
                {
                    if(norm2(g) < norm2(g_old))
                    {
                        x = x_trial; 
                        taken = 1; 
                        g_old = g; 
                        inner_counter = 0;
                    }
                    else
                    {
                        taken=0;
                        inner_counter++; 
                    }


                    // clamping values of tau
                    if(std::isinf(tau))
                        tau = 9e249;        
                    else if(std::isnan(tau))
                        tau = 1e-10;
                    else if(tau==0)
                        tau = 1e-10;                                        
                }

                r_norm = dot(s,g); 
                g_norm = norm2(g);

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm, tau, nu,  double(taken)});

                r_norm = 9e9; 
                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, r_norm, s_norm);
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
