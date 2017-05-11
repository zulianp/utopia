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
        typedef UTOPIA_SCALAR(Vector)    Scalar;

        public:
            PETSCUtopiaNonlinearFunction(SNES snes, const Vector & x_init = local_zeros(1), const Vector & rhs = local_zeros(1)) 
            :
                Function<Matrix, Vector>(x_init, rhs),
                snes_(snes), 
                first_grad(0), 
                first_energy(0)
            {
            
            }


            virtual bool gradient(const Vector &x, Vector &g) const override
            {
                // initialization of gradient vector... 
                if(empty(g)){
                    g  = local_zeros(local_size(x));; 
                }
                
                SNESComputeFunction(snes_, raw_type(x), raw_type(g));   

                // if(local_size(g)==local_size(this->_rhs))
                // {
                //    std::cout<<"grad:: yes rhs ... \n"; 
                //     g = g - this->_rhs; 
                // }

                return true; 
            }

            virtual bool hessian(const Vector &x, Matrix &hessian) const override
            {
                // std::cout<<"before hessian \n"; 
                SNESComputeJacobian(snes_, raw_type(x), snes_->jacobian,  snes_->jacobian_pre);
                //hessian = sparse_mref(snes_->jacobian); 
                convert(snes_->jacobian, hessian); 
                // std::cout<<"after hessian \n"; 

                return true; 
            }

            virtual bool value(const Vector &x, typename Vector::Scalar &result) const override 
            {

                // if(first_energy == 0) 
                // {
                //     result = 9e12; 
                //     first_energy++; 
                //     return true; 
                // }

                // hack to have fresh energy (MOOSE post-processor does things in strange way )
                Vector grad = 0 * x; 
                this->gradient(x, grad);         

            

                SNESComputeObjective(snes_, raw_type(x), &result); 

               // result = 0.5 * norm2(grad) * norm2(grad); 



                // Matrix H;
                // this->hessian(x, H);    

                // Scalar EN2 = 0.5 * dot(x, H * x);
                // Scalar EN3 = EN2 - dot(x, grad);

               //  std::cout<<"multiplication result: "<< EN2 << "      E3:    "<< EN3 <<  "      result:    "<< result << "   \n"; 



                // result -= dot(x,grad); 


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