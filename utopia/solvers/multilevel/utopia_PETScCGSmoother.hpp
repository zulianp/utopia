/*
* @Author: alenakopanicakova
* @Date:   2016-06-27
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-11
*/

#ifndef UTOPIA_PETSC_CG_SMOOTHER_HPP
#define UTOPIA_PETSC_CG_SMOOTHER_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"

// extern "C" 
// {
#include "petscmat.h"
#include "petscvec.h"
#include <petscksp.h>
// }


    namespace utopia 
    {
        /**
         * @brief     Wrapper for PETSc CG to be used as smoother. 
         *  
         * @tparam     Matrix  
         * @tparam     Vector  
         */
        template<class Matrix, class Vector>
        class Petsc_CG_Smoother : public  Smoother<Matrix, Vector>,  public IterativeSolver<Matrix, Vector>
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
            typedef utopia::IterativeSolver<Matrix, Vector> Solver;
            typedef utopia::Smoother<Matrix, Vector> Smoother;

        public:
        Petsc_CG_Smoother(const Parameters params = Parameters()) 
        {
            set_parameters(params); 
        }

        // /**
        //  *
        //  *   CG smoother smoothing
        //  *
        //  */
        bool smooth(const Matrix &A, const Vector &rhs, Vector &x) override
        {        
            KSP solver;
            KSPCreate(A.implementation().communicator(), &solver);
            KSPSetOperators(solver, raw_type(A), raw_type(A));
            
            // KSPSetType(solver, KSPCGS);
            KSPSetType(solver, KSPBCGS);
            KSPSetInitialGuessNonzero(solver, PETSC_TRUE);

            KSPSetUp(solver);
            KSPSetTolerances(solver, this->rtol(), this->atol(), this->stol(), this->sweeps());
            KSPSolve(solver, raw_type(rhs), raw_type(x));

            return true; 
        }

        bool solve(const Matrix &A, const Vector &rhs, Vector &x) override
        {
            this->init_solver("Petsc CG", {"it. ", "||r||" }); 

            KSP solver;
            KSPCreate(A.implementation().communicator(), &solver);
            KSPSetOperators(solver, raw_type(A), raw_type(A));
            
            // KSPSetType(solver, KSPCGS);
            KSPSetType(solver, KSPBCGS);
            KSPSetInitialGuessNonzero(solver, PETSC_TRUE);

            KSPSetUp(solver);
            KSPSetTolerances(solver, this->rtol(), this->atol(), this->stol(),  this->max_it());
            KSPSolve(solver, raw_type(rhs), raw_type(x));

            return true; 
        }


        void set_parameters(const Parameters params) override
        {
            Smoother::set_parameters(params); 
            Solver::set_parameters(params); 
        }


};

}

#endif //UTOPIA_PETSC_CG_SMOOTHER_HPP

