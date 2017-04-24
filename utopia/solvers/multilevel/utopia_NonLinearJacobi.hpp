/*
* @Author: alenakopanicakova
* @Date:   2017-04-17
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-04-22
*/

#ifndef UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP
#define UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"

#include "utopia_PETScFunction.hpp"



namespace utopia 
{

    /**
     * @brief      Nonlinear Jacobi smoother. Used in nonlinear MG as well as in FAS.
     *  
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class NonLinearJacobi: public NonLinearSmoother<Matrix, Vector> 
    {
            typedef UTOPIA_SCALAR(Vector)                           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
            typedef utopia::NonLinearSmoother<Matrix, Vector>       Smoother;
            typedef utopia::Function<Matrix, Vector>                Function;

        public:
        NonLinearJacobi(const Parameters params = Parameters()) 
        { 
            set_parameters(params); 
        }

        virtual void set_parameters(const Parameters params) override
        {
            Smoother::set_parameters(params); 
        }



        // takes function
        // - take PETSC nonlinear function - to get SNES ??? 
        // dynamic cast ?? 
        virtual bool nonlinear_smooth(Function & fun,  Vector &x, const Vector &rhs) override
        {
            Vector F = local_zeros(local_size(x));
               
            for(int i =0; i < this->sweeps(); i++)
            {
                Matrix D; 
                fun.hessian(x, D);                 
                fun.gradient(x, F); 
         
                Vector d = 1/diag(D); 
                D = diag(d); 
         
                x = x - this->relaxation_parameter() * (D * (F - rhs)); 
            }
            
            return true; 

        }
    };

}

#endif //UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP

