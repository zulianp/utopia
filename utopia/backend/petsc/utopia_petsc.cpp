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
#include "utopia_ProjectedGaussSeidel_impl.hpp"
#include "utopia_RMTR.hpp"
#include "utopia_RMTRVcycleImpl.hpp"
#include "utopia_RMTR_inf.hpp"
#include "utopia_SaddlePoint.hpp"
#include "utopia_SemismoothNewton_impl.hpp"
#include "utopia_TRBoundsGratton.hpp"
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_impl.hpp"

// explicit instantiations
namespace utopia {

    template class ConjugateGradient<PetscMatrix, PetscVector>;
    template class ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE>;
    template class GaussSeidel<PetscMatrix, PetscVector, PETSC>;
    template class SPBlockConjugateGradient<PetscMatrix, PetscVector>;
    template class BiCGStab<PetscMatrix, PetscVector, HOMEMADE>;

    // qp-solvers
    template class ProjectedGaussSeidel<PetscMatrix, PetscVector>;
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

    void optimize_nnz(PetscMatrix &A) {
        auto rr = row_range(A);
        auto cr = col_range(A);
        auto ls = local_size(A);
        auto gs = size(A);

        std::vector<PetscInt> d_nnz(rr.extent(), 0), o_nnz(rr.extent(), 0);
        A.read([&](const utopia::SizeType i, const utopia::SizeType j, const PetscScalar val) {
            if (std::abs(val) > 1e-18) {
                if (cr.inside(j)) {
                    ++d_nnz[i - rr.begin()];
                } else {
                    ++o_nnz[i - rr.begin()];
                }
            }
        });

        PetscMatrix A_opt;

        A_opt.matij_init(A.communicator(), A.type_override(), ls.get(0), ls.get(1), gs.get(0), gs.get(1), d_nnz, o_nnz);

        {
            Write<PetscMatrix> w_A(A_opt);
            A.read([&](const SizeType i, const SizeType j, const PetscScalar val) {
                if (std::abs(val) > 1e-18) {
                    A_opt.set(i, j, val);
                }
            });
        }

        A = std::move(A_opt);
    }

    bool is_diagonally_dominant(const PetscMatrix &A) {
        PetscVector d = diag(A);
        PetscVector o(layout(d));

        {
            Write<PetscVector> w_o(o);
            A.read([&o](const SizeType i, const SizeType j, const PetscScalar val) {
                if (i != j) {
                    o.add(i, std::abs(val));
                }
            });
        }

        PetscVector diff = d - o;
        PetscScalar m = min(diff);
        return m > 0.;
    }

    void local_block_view(const PetscMatrix &mat, PetscMatrix &block) {
        Mat M;
        auto ierr = MatGetDiagonalBlock(mat.raw_type(), &M);
        assert(ierr == 0);
        UTOPIA_UNUSED(ierr);

        block.wrap(M);
        block.update_mirror();
    }
}  // namespace utopia
