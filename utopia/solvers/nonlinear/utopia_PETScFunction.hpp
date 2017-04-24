// /*! \file passo_Newton_Solver.cpp
//     Passo nonlinear function
//     in order to have interface between petsc and utopia
//     Created by Alena Kopanicakova 
// */

#ifndef PETSC_BASED_UTOPIA_NONLINEAR_FUNCTION_HPP
#define PETSC_BASED_UTOPIA_NONLINEAR_FUNCTION_HPP


#include <petscsnes.h>
#include <petsc/private/snesimpl.h>
#include "utopia_Core.hpp"
#include "utopia_PETScTypes.hpp"

namespace utopia 
{
    template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class PETSCUtopiaNonlinearFunction {};
      
    
    template<class Matrix, class Vector>
    class PETSCUtopiaNonlinearFunction<Matrix, Vector, PETSC> : public Function<Matrix, Vector> 
    {

        public:
            PETSCUtopiaNonlinearFunction(SNES snes) 
            : snes_(snes), first_grad(0), first_energy(0)
            {
            
            }



            virtual bool gradient(const Vector &x, Vector &g) const override
            {
                // initialization of gradient vector... 
                if(empty(g)){
                    g  = local_zeros(local_size(x));; 
                }
                
                SNESComputeFunction(snes_, raw_type(x), raw_type(g));   

                if(local_size(g)==local_size(this->_rhs))
                {
                    g = g - this->_rhs; 
                }

                return true; 
            }

            virtual bool hessian(const Vector &x, Matrix &hessian) const override
            {

                SNESComputeJacobian(snes_, raw_type(x), snes_->jacobian,  snes_->jacobian_pre);
                hessian = sparse_mref(snes_->jacobian); 

                return true; 
            }

            virtual bool value(const Vector &x, typename Vector::Scalar &result) const override 
            {

                // // hack to have fresh energy (MOOSE post-processor does things in strange way )
                // Vector grad = 0 * x; 
                // this->gradient(x, grad); 

                if(first_energy == 0) 
                {
                    result = 9e12; 
                    first_energy++; 
                    return true; 
                }

                SNESComputeObjective(snes_, raw_type(x), &result); 
                return true; 
            }

            virtual void getSNES(SNES &snes)
            {
                snes = snes_; 
            }


        private:

            SNES snes_;
            mutable PetscInt first_grad; 
            mutable PetscInt first_energy; 
        };

    }



#endif  //PETSC_BASED_UTOPIA_NONLINEAR_FUNCTION_HPP