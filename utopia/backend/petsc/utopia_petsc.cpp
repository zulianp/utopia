// https://fossies.org/diffs/petsc/3.7.7_vs_3.8.0/
// https://www.mcs.anl.gov/petsc/documentation/changes/39.html

#include "utopia_petsc.hpp"
#include "utopia_AffineSimilarity.hpp"
#include "utopia_BiCGStab_impl.hpp"
#include "utopia_BlockQPSolver_impl.hpp"
#include "utopia_ConjugateGradient_impl.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FAS.hpp"
#include "utopia_MG_OPT.hpp"
#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_RMTR.hpp"
#include "utopia_RMTRVcycleImpl.hpp"
#include "utopia_RMTR_inf.hpp"
#include "utopia_SaddlePoint.hpp"
#include "utopia_SemismoothNewton_impl.hpp"
#include "utopia_TRBoundsGratton.hpp"
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_impl.hpp"

#include "utopia_ILU_impl.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"
#include "utopia_petsc_ILUDecompose.hpp"

// explicit instantiations
namespace utopia {

    template class ConjugateGradient<PetscMatrix, PetscVector>;
    template class ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE>;
    template class GaussSeidel<PetscMatrix, PetscVector, PETSC>;
    template class SPBlockConjugateGradient<PetscMatrix, PetscVector>;
    template class BiCGStab<PetscMatrix, PetscVector, HOMEMADE>;
    template class ILU<PetscMatrix, PetscVector>;

    // qp-solvers
    template class BlockQPSolver<PetscMatrix, PetscVector>;
    template class SemismoothNewton<PetscMatrix, PetscVector, PETSC>;
    template class SemismoothNewton<PetscMatrix, PetscVector, HOMEMADE>;

    // petsc non-linear solvers
    template class NonLinearGaussSeidel<PetscMatrix, PetscVector>;
    template class Multigrid<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>;

    template class RMTR_l2<PetscMatrix, PetscVector, FIRST_ORDER>;
    template class RMTR_inf<PetscMatrix, PetscVector, TRBoundsGratton<PetscMatrix, PetscVector> >;

    template class FAS<PetscMatrix, PetscVector>;
    template class MG_OPT<PetscMatrix, PetscVector>;

    template class AffineSimilarity<PetscMatrix, PetscVector>;

}  // namespace utopia
