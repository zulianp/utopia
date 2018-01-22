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
#include "utopia_petsc_KSPSolver.hpp"

// extern "C" 
// {
#include "petscmat.h"
#include "petscvec.h"
#include <petscksp.h>
// }

namespace utopia  {
    /**
     * @brief     Wrapper for Petsc CG to be used both as smoother and solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    
    template<typename Matrix, typename Vector>
    class ConjugateGradient<Matrix, Vector, PETSC> : public Smoother<Matrix, Vector>, public KSPSolver<Matrix, Vector, PETSC> {
    public:
        ConjugateGradient(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : KSPSolver<Matrix, Vector, PETSC>(params), preconditioner_(preconditioner)
        {
            this->pc_type(preconditioner_);
            this->ksp_type("cg");
        }
        
        inline void set_parameters(const Parameters params) override
        {
            Parameters params_copy = params;
            params_copy.lin_solver_type("cg");
            params_copy.preconditioner_type(preconditioner_.c_str());
            KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
            Smoother<Matrix, Vector>::set_parameters(params);
        }
        
        inline void pc_type(const std::string & preconditioner) override
        {
            preconditioner_ = preconditioner;
            KSPSolver<Matrix, Vector, PETSC>::pc_type(preconditioner_);
        }
        
        //FIXME use parameters from KSPSolver
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
            
            KSPDestroy(&solver);
            return true;
        }
        
    private:
        std::string preconditioner_;
    };
    
}

#endif //UTOPIA_PETSC_CG_SMOOTHER_HPP

