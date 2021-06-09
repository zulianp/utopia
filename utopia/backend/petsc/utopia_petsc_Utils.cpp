#include "utopia_petsc_Utils.hpp"

#include "utopia_petsc_Matrix_impl.hpp"

#include "utopia_Diag.hpp"
#include "utopia_Writable.hpp"

#include "utopia_Core.hpp"

#include "utopia_RowView.hpp"
#include "utopia_petsc_RowView.hpp"

namespace utopia {

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

    // void local_view(const PetscVector &vec, PetscVector &lv) {
    //     assert(false);
    //     UTOPIA_UNUSED(vec);
    //     UTOPIA_UNUSED(lv);
    //     // VecGetLocalVector(Vec v,Vec w)
    // }

}  // namespace utopia
