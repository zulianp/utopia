#ifndef UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP
#define UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"



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

        virtual bool nonlinear_smooth(Function & fun,  Vector &x, const Vector &rhs) override
        {
            Vector g = local_zeros(local_size(x));
               
            for(int i =0; i < this->sweeps(); i++)
            {
                Matrix H; 
                fun.hessian(x, H);             
                fun.gradient(x, g); 
                g -= rhs; 
         
                Vector d = 1./diag(H); 
                H = diag(d); 

                x = x - (this->damping_parameter() * H * g); 
            }
            
            return true; 

        }
    };

}

#endif //UTOPIA_NONLINEAR_JACOBI_SMOOTHER_HPP

