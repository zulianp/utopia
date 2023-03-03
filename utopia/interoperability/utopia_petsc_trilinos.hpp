#ifndef UTOPIA_PETSC_TRILINOS_HPP
#define UTOPIA_PETSC_TRILINOS_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#ifdef UTOPIA_ENABLE_PETSC

#include "utopia_CrossBackendLinearSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_Factorization.hpp"
#include "utopia_petsc_GaussSeidel.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_KSPSolvers.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_trilinos_Types.hpp"

#include "utopia_petsc_trilinos_ConvertTensor.hpp"

namespace utopia {

    template <typename Matrix, typename Vector>
    class KSPSolver<Matrix, Vector, TRILINOS>
        : public CrossBackendLinearSolverAndSmoother<Matrix,
                                                     Vector,
                                                     PetscMatrix,
                                                     PetscVector,
                                                     KSPSolver<PetscMatrix, PetscVector, PETSC> > {};

    // FIXME remove me once the belos solver works
    template <typename Matrix, typename Vector>
    class Factorization<Matrix, Vector, TRILINOS>
        : public CrossBackendLinearSolver<Matrix,
                                          Vector,
                                          PetscMatrix,
                                          PetscVector,
                                          Factorization<PetscMatrix, PetscVector, PETSC> > {};

    // FIXME remove me once the belos solver works
    template <typename Matrix, typename Vector>
    class BiCGStab<Matrix, Vector, TRILINOS>
        : public CrossBackendLinearSolverAndSmoother<Matrix,
                                                     Vector,
                                                     PetscMatrix,
                                                     PetscVector,
                                                     BiCGStab<PetscMatrix, PetscVector, PETSC> > {};

    // FIXME remove me once the belos solver works
    template <typename Matrix, typename Vector>
    class GaussSeidel<Matrix, Vector, TRILINOS>
        : public CrossBackendLinearSolverAndSmoother<Matrix,
                                                     Vector,
                                                     PetscMatrix,
                                                     PetscVector,
                                                     GaussSeidel<PetscMatrix, PetscVector, PETSC> > {};

    // FIXME remove me once the belos solver works
    template <typename Matrix, typename Vector>
    class GMRES<Matrix, Vector, TRILINOS>
        : public CrossBackendLinearSolverAndSmoother<Matrix,
                                                     Vector,
                                                     PetscMatrix,
                                                     PetscVector,
                                                     GMRES<PetscMatrix, PetscVector, PETSC> > {};

}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS
#endif  // UTOPIA_ENABLE_PETSC

#endif  // UTOPIA_PETSC_TRILINOS_HPP
