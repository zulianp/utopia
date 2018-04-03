/*
* @Author: alenakopanicakova
* @Date:   2017-04-17
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/

#ifndef UTOPIA_NONLINEAR_PETSC_GS_HPP
#define UTOPIA_NONLINEAR_PETSC_GS_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"

#include "utopia_petsc_SNESFunction.hpp"

#include <petsc/private/snesimpl.h>
#include "petscsnes.h"  



namespace utopia {
 //FIXME move it to utopia/solvers/smoothers
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

            std::cerr<<"------- utopia_petsc_NonLinearGaussSeidl::nonlinear_smooth:: there is something wrong int this function -------------- \n"; 
            // TODO:: understand problem on snes side... 


            if(dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun))
            {
                PETSCUtopiaNonlinearFunction<Matrix, Vector> * fun_petsc = dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun);
                
                SNES snes;
                fun_petsc->getSNES(snes); 

                SNESSetFromOptions(snes); 

                // weird .......
                // SNESSetType(snes, SNESNRICHARDSON);
                

                SNESType type; 
                SNESGetType(snes, &type); 

                std::cout<<"snes - type : "<< type << "      --- \n"; 

                Vector y = x;  

                SNESComputeJacobian(snes, raw_type(x), snes->jacobian,  snes->jacobian_pre);
                SNESComputeFunction(snes, raw_type(x), raw_type(y)); 


                // SNES pc; 
                // if(!this->verbose())
                //     SNESMonitorCancel(snes);
                                    
                                // SNESNRICHARDSON
                                // SNESNGS
                                // SNESNGMRES 
                                // SNESNCG
                SNESSetType(snes, SNESNGS );
                SNESSetTolerances(snes, 0.0, 0.0, 0.0, this->sweeps(), PETSC_DEFAULT);

                SNESLineSearch linesearch; 
                SNESGetLineSearch(snes, &linesearch);
                SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC); 

                SNESSolve(snes, raw_type(rhs), raw_type(x)); 
                snes->vec_rhs =  NULL; 
                
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

