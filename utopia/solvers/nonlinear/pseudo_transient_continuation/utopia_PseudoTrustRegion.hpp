#ifndef UTOPIA_PSEUDO_TRUST_REGION_HPP
#define UTOPIA_PSEUDO_TRUST_REGION_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_EigenSolver.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    
    template<class Matrix, class Vector>
    class PseudoTrustRegion : public NewtonBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
        typedef utopia::EigenSolver<Matrix, Vector>         EigenSolver;

        using NewtonBase<Matrix, Vector>::print_statistics; 


    public:
       PseudoTrustRegion(   const std::shared_ptr <Solver> &linear_solver,
                            const std::shared_ptr<EigenSolver> & eigen_solver,  
                            const Parameters params                       = Parameters() ):
                            NewtonBase<Matrix, Vector>(linear_solver, params), 
                            eigen_solver_(eigen_solver), 
                            eps_(1e-12)
                            {

                            }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

            if(this->verbose())
                this->init_solver("PseudoTrustRegion", {" it. ", "|| g ||", " E ",  "rho", "tau", "|| Delta x || "});

            bool converged = false; 
            SizeType it = 0; 

            Scalar g_norm, s_norm=9e9, tau, lambda_min;
            Scalar ared, pred, rho; 

            Vector g = local_zeros(local_size(x)), eigenvector_min, s, x_trial; 
            Matrix H, H_dumped; 

            Matrix I = local_identity(local_size(x).get(0), local_size(x).get(0)); 

            fun.gradient(x, g); 
            g_norm = norm2(g); 

            Scalar energy_old, energy_new, energy; 
            fun.value(x, energy_new);

            tau = 1; 
        
            while(!converged)
            {
                fun.hessian(x, H); 
                fun.gradient(x, g); 

                H_dumped = H + 1./tau * I; 

                eigen_solver_->portion_of_spectrum("smallest_real"); 
                eigen_solver_->number_of_eigenvalues(1); 
                eigen_solver_->solve(H_dumped); 
                eigen_solver_->get_real_eigenpair(0, lambda_min, eigenvector_min); 

                if(lambda_min >= eps_)
                {
                    s = 0 * x; 
                    this->linear_solve(H_dumped, -1.0 * g, s);
                    x_trial = x + s; 

                    fun.value(x, energy_old); 
                    fun.value(x_trial, energy_new); 
                    ared = energy_old - energy_new; 

                    pred = (-1.0 * dot(g, s) -0.5 *dot(H* s, s));
                    rho = ared/pred; 

                    if(rho < 1./4.)
                        tau *=0.5; 
                    else if(rho > 3./4.)
                        tau *=2.0;
                }
                else
                {
                    rho = -1.0; 
                    tau  *= 0.5; 
                }


                if(rho > 0.0)
                {
                    x = x_trial; 
                }

                fun.gradient(x, g); 
                g_norm = norm2(g); 
                s_norm = norm2(s); 
                fun.value(x, energy); 

                it++; 

                if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm, energy,  rho, tau, s_norm});

                converged = this->check_convergence(it, g_norm, 9e9, s_norm);


            } // outer solve loop while(!converged)


            return true;
        }
    
    private: 
        std::shared_ptr<EigenSolver> eigen_solver_;     
        Scalar eps_; 


    };

}
#endif //UTOPIA_PSEUDO_TRUST_REGION_HPP
