#ifndef UTOPIA_PETSC_HPP
#define UTOPIA_PETSC_HPP

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Solver_Traits.hpp"
#include "utopia_petsc_Traits.hpp"

#include "utopia_petsc_Error.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_petsc_Types.hpp"

#include "utopia_petsc_RowView.hpp"

#include "utopia_petsc_LinearSolverFactory.hpp"
#include "utopia_petsc_TrustRegionFactory.hpp"
#include "utopia_petsc_solvers.hpp"

//#include "utopia_petsc_Each.hpp"
#include "utopia_petsc_Newton.hpp"
#include "utopia_petsc_TaoTRQP.hpp"

#ifdef UTOPIA_WITH_SLEPC
#include "utopia_petsc_Slepc.hpp"
#endif

// very much experimental files for the moment
//#include "utopia_petsc_Each.hpp"
#include "utopia_petsc_SNES.hpp"
#include "utopia_petsc_build_ksp.hpp"
#include "utopia_petsc_debug.hpp"

#include "utopia_petsc_Eval.hpp"
#include "utopia_petsc_Layout.hpp"

/// FIXME

// #include "utopia_petsc_Matrix_impl.hpp"
// #include "utopia_petsc_Vector_impl.hpp"

namespace utopia {
    void optimize_nnz(PetscMatrix &A);
    bool is_diagonally_dominant(const PetscMatrix &A);
    void local_block_view(const PetscMatrix &mat, PetscMatrix &block);
}  // namespace utopia

#endif  // UTOPIA_PETSC_HPP
