/*
* @Author: alenakopanicakova
* @Date:   2016-07-14
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-07-29
*/

#ifndef PETSC_UTOPIA_SNES_NONLINEAR_FUNCTION_HPP
#define PETSC_UTOPIA_SNES_NONLINEAR_FUNCTION_HPP


#include <petscsnes.h>
#include <petsc/private/snesimpl.h>
#include <utopia.hpp>

namespace utopia 
{
    /**
     * @brief      Function carries over SNES application context for nonlinear solve \n.
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector> 
    class PETSCUtopiaNonlinearFunction : public Function<Matrix, Vector> 
    {

        public:
            PETSCUtopiaNonlinearFunction(SNES snes) 
            : snes(snes), first_grad(0) 
            {
            
            }


            bool gradient(const Vector &x, Vector &gradient) const override
            {
                if(first_grad == 0) 
                {
                    utopia::convert(snes->vec_sol_update, gradient);
                    first_grad++;
                }

                SNESComputeFunction(snes,utopia::raw_type(x), utopia::raw_type(gradient));   
                return true; 
            }

            bool hessian(const Vector &x, Matrix &hessian) const override
            {
                SNESComputeJacobian(snes, utopia::raw_type(x), snes->jacobian,  snes->jacobian_pre);
                utopia::convert(snes->jacobian, hessian);
                return true; 
            }

            bool value(const Vector &point, typename Vector::Scalar &result) const override 
            {

                // merit ...
                SNESComputeFunction(snes,utopia::raw_type(point), snes->vec_sol_update);   
                PetscReal val; 
                VecNorm(snes->vec_sol_update, NORM_2, &val);
                result = val * val; 

                // energy ...
                // SNESComputeObjective(snes, utopia::raw_type(point), &val); 
                // result = val; 
                return true; 
            }

            bool residual(const Vector &point, Vector &residual)  override 
            {

                if(first_resid == 0) 
                {
                    utopia::convert(snes->vec_sol_update, residual);
                    first_resid++;
                }

                SNESComputeFunction(snes,utopia::raw_type(point), utopia::raw_type(residual));  

                return true; 
            }


            bool jacobian(const Vector &x, Matrix &jacobian) override
            {

                SNESComputeJacobian(snes, utopia::raw_type(x), snes->jacobian,  snes->jacobian_pre);
                utopia::convert(snes->jacobian_pre, jacobian);
                return true; 
            }


        private:

            SNES snes;
            mutable PetscInt first_grad = 0 ;
            mutable PetscInt first_resid = 0;
        };

    }



#endif  //PETSC_UTOPIA_SNES_NONLINEAR_FUNCTION_HPP