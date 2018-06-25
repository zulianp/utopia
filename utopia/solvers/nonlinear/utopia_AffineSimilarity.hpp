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

            Vector g, s, rhs, x_trial, x0, g0;
            Matrix H, A;

            Scalar g_norm, g0_norm=9e9, r_norm=9e9, s_norm=9e9, tau;
            SizeType it = 0;

            bool converged = false, taken=0;

            x0 = x; 

            fun.gradient(x, g);
            g0 = -1.0*g;

            g0_norm = norm2(g0);
            g_norm = g0_norm;


            Scalar L0 = g0_norm; 

            this->init_solver("Affine similarity", {" it. ", "|| g ||", "<s,g>", "|| s_k || ", "tau", "nu", "nu_test",   "taken"});

            if(this->verbose_)
                PrintInfo::print_iter_status(it, {g_norm, 1, 0});
            it++;

            tau = norm2(g); 
            tau = 1/tau; 


            fun.hessian(x, H);
            Matrix H0 = -1.0* H; 

            Scalar nu_zero_approx = dot(g0, H0*g0)/(g0_norm*g0_norm); 

            std::cout<<"nu_zero_approx: "<< nu_zero_approx << "  \n"; 

            Vector g_old = g0; 

            while(!converged)
            {
                //find direction step
                s = local_zeros(local_size(x));
                fun.hessian(x, H);
                H *= -1.0; 

                fun.gradient(x, g);                
                g *= -1.0;                 
                
                A = (M_ - tau*H); 
                rhs = tau*g; 
                

                this->linear_solve(A, rhs, s);

                x_trial = x+s; 

                fun.gradient(x_trial, g);  
                g *= -1.0;     

                s = 1./tau * s; 
                Vector gs_diff = g-s; 
                s_norm = norm2(s); 

                Scalar nu = dot(s, s-g_old)/(s_norm*s_norm);
                Scalar nu_test = tau * dot(s, H0*s)/ (s_norm*s_norm); 

                {
                    if( norm2(g) < g0_norm)
                    {
                        x = x_trial; 
                        taken = 1; 
                    }

                    
                    Scalar L2 = 2.0 * norm2(gs_diff); 
                    Scalar help = (tau*tau * s_norm*s_norm); 
                    L2 = L2/ help;
                    
                    Scalar help2 = L2 * s_norm; 
                    tau = std::abs(nu)/help2; 

                    // clamping taus
                    if(std::isinf(tau))
                        tau = 9e249;                     
                }

                r_norm = dot(s,g); 
                g_old = g; 

                fun.gradient(x, g);     
                g *= -1.0;          
                g_norm = norm2(g);


                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm, tau, nu, nu_test,  double(taken)});


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
