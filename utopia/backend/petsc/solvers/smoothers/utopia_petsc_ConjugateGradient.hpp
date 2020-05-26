#ifndef UTOPIA_PETSC_CG_SMOOTHER_HPP
#define UTOPIA_PETSC_CG_SMOOTHER_HPP

#include "utopia_ConjugateGradient.hpp"
#include "utopia_Core.hpp"
#include "utopia_petsc_KSPSolver.hpp"

#include <petscksp.h>
#include "petscmat.h"
#include "petscvec.h"

namespace utopia {
    /**
     * @brief     Wrapper for Petsc CG to be used both as smoother and solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <typename Matrix, typename Vector>
    class ConjugateGradient<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        ConjugateGradient(const std::string &preconditioner = "jacobi") : KSPSolver<Matrix, Vector, PETSC>() {
            this->pc_type(preconditioner);
            this->ksp_type("cg");
        }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_CG_SMOOTHER_HPP
