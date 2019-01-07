#ifndef UTOPIA_PETSC_CG_SMOOTHER_HPP
#define UTOPIA_PETSC_CG_SMOOTHER_HPP

#include "utopia_Core.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_ConjugateGradient.hpp"

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
        ConjugateGradient(const std::string &preconditioner = "jacobi")
        : KSPSolver<Matrix, Vector, PETSC>()
        {
            this->pc_type(preconditioner);
            this->ksp_type("cg");
        }
    };
}

#endif //UTOPIA_PETSC_CG_SMOOTHER_HPP

