/*
* @Author: alenakopanicakova
* @Date:   2017-04-17
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-13
*/

#ifndef UTOPIA_NONLINEAR_PETSC_GS_HPP
#define UTOPIA_NONLINEAR_PETSC_GS_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"


#ifdef WITH_PETSC
    #include "utopia_PETScFunction.hpp"
    #include <petsc/private/snesimpl.h>
    #include "petscsnes.h"  
#endif //WITH_PETSC


namespace utopia 
{
  template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  class NonLinearGaussSeidel {};


        /**
         * @brief      Wrapper for PETSC implementation of NGS. 
         * @tparam     Matrix  
         * @tparam     Vector  
         */
    template<class Matrix, class Vector>
    class NonLinearGaussSeidel<Matrix, Vector, PETSC> : public NonLinearSmoother<Matrix, Vector> 
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
        typedef utopia::NonLinearSmoother<Matrix, Vector>       Smoother;
        typedef utopia::Function<Matrix, Vector>                Function;

        public:
        NonLinearGaussSeidel(const Parameters params = Parameters()) 
        { 
            set_parameters(params); 
        }

        virtual void set_parameters(const Parameters params) override
        {
            Smoother::set_parameters(params); 
        }

        virtual bool nonlinear_smooth(Function & fun,  Vector &x, const Vector &rhs) override
        {

            if(dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun))
            {
                PETSCUtopiaNonlinearFunction<Matrix, Vector> * fun_petsc = dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun);
                
                SNES snes; 
                fun_petsc->getSNES(snes); 

                PetscScalar result; 
                // SNESComputeObjective(snes, raw_type(x), &result); 
                // std::cout<<"before GS smoother, E:  "<< result <<  "  \n"; 

                SNESSetFromOptions(snes); 
                SNESSetType(snes, SNESNRICHARDSON);

                SNESComputeJacobian(snes, raw_type(x), snes->jacobian,  snes->jacobian_pre);

                // SNES pc; 
                SNESSetType(snes, SNESNGS);
                SNESSetTolerances(snes, 0.0, 0.0, 0.0, this->sweeps(), PETSC_DEFAULT);

                SNESLineSearch linesearch; 
                SNESGetLineSearch(snes, &linesearch);
                SNESLineSearchSetType(linesearch, SNESLINESEARCHL2); 

                SNESSolve(snes, raw_type(rhs), raw_type(x)); 
                snes->vec_rhs =  NULL; 

                // SNESComputeObjective(snes, raw_type(x), &result); 
                // std::cout<<"after GS smoother, E:  "<< result <<  "  \n"; 
                
            }
            else
            {
                std::cout<<"utopia:: PETSC NONLINEAR GS works just with PETSC function.... \n"; 
                exit(1);
            }
            return true; 
        }
    };

}

#endif //UTOPIA_NONLINEAR_PETSC_GS_HPP

