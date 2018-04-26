#ifndef UTOPIA_PETSC_CG_SMOOTHER_HPP
#define UTOPIA_PETSC_CG_SMOOTHER_HPP

#include "utopia_Core.hpp"
#include "utopia_petsc_KSPSolver.hpp"

#include "petscmat.h"
#include "petscvec.h"
#include <petscksp.h>

namespace utopia  {
    /**
     * @brief     Wrapper for Petsc CG to be used both as smoother and solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<typename Matrix, typename Vector>
    class ConjugateGradient<Matrix, Vector, PETSC> : /*public Smoother<Matrix, Vector>,*/ public KSPSolver<Matrix, Vector, PETSC> {
    public:
        ConjugateGradient(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : KSPSolver<Matrix, Vector, PETSC>(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("cg");
        }
        
        inline void set_parameters(const Parameters params) override
        {
            Parameters params_copy = params;
            params_copy.lin_solver_type("cg");
            params_copy.preconditioner_type(this->pc_type().c_str());
            KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
            Smoother<Matrix, Vector>::set_parameters(params);
        }
    };
}

#endif //UTOPIA_PETSC_CG_SMOOTHER_HPP

