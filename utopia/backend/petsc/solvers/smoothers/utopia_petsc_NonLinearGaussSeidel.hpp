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
#include "utopia_petsc_SNES.hpp"

#include <petsc/private/snesimpl.h>
#include "petscsnes.h"  



namespace utopia {
 //FIXME move it to utopia/solvers/smoothers
  template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  class NonLinearGaussSeidel {};




    template<class Matrix, class Vector>
    class NonLinearGaussSeidel<Matrix, Vector, PETSC> : public SNESSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

        typedef utopia::SNESSolver<Matrix, Vector>                  SNESSolver;
        typedef utopia::Function<Matrix, Vector>                    Function;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver    LinearSolver;


        public:
        NonLinearGaussSeidel(   const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(), 
                                const Parameters params = Parameters()) 
                                : SNESSolver(linear_solver, params)
        { 
            set_parameters(params); 
            this->set_snes_type("ngs"); 

        }

        virtual void set_parameters(const Parameters params) override
        {
            SNESSolver::set_parameters(params); 
        }

    protected: 
        virtual void set_snes_options(SNES & snes,  const Scalar & atol     = SNESSolver::atol(), 
                                                    const Scalar & rtol     = SNESSolver::rtol(), 
                                                    const Scalar & stol     = SNESSolver::stol(), 
                                                    const SizeType & max_it = SNESSolver::max_it()) override 
        {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it); 

            // we need to allocate hessian for coloring computation
            SNESComputeJacobian(snes, snes->vec_sol, snes->jacobian,  snes->jacobian_pre);

            SNESLineSearch linesearch; 
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC); 
        }

    };

}

#endif //UTOPIA_NONLINEAR_PETSC_GS_HPP

