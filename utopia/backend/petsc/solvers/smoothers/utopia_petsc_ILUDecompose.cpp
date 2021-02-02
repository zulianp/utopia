#include "utopia_petsc_ILUDecompose.hpp"
#include "utopia_petsc.hpp"

#include <vector>
#include "utopia_Views.hpp"

#include "utopia_CRSToBlockCRS.hpp"
#include "utopia_ILUDecompose.hpp"

namespace utopia {

    class PetscSeqAIJRaw {
    public:
        PetscSeqAIJRaw(Mat raw_mat) : raw_mat(raw_mat) {
            err = MatGetRowIJ(raw_mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
            assert(err == 0);
            assert(done == PETSC_TRUE);

            if (!done) {
                utopia::err()
                    << "PetscMatrix::read_petsc_seqaij_impl(const Op &op): MatGetRowIJ failed to provide what "
                       "was asked.\n";
                abort();
            }

            MatSeqAIJGetArray(raw_mat, &array);
        }

        ~PetscSeqAIJRaw() {
            MatSeqAIJRestoreArray(raw_mat, &array);
            err = MatRestoreRowIJ(raw_mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
            assert(err == 0);
            UTOPIA_UNUSED(err);
        }

        Mat raw_mat;
        PetscInt n = 0;
        const PetscInt *ia{nullptr};
        const PetscInt *ja{nullptr};
        PetscScalar *array{nullptr};
        PetscBool done;
        PetscErrorCode err{0};
    };

    void ILUDecompose<PetscMatrix, PETSC>::decompose(const PetscMatrix &mat, PetscMatrix &out, const bool modified) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscMatrix l_mat;
        local_block_view(mat, l_mat);

        // perform copy
        out.copy(l_mat);

        PetscSeqAIJRaw m_raw(out.raw_type());
        PetscInt n = m_raw.n;
        PetscInt nnz = m_raw.ia[n];

        const PetscInt *ia = m_raw.ia;
        const PetscInt *ja = m_raw.ja;
        PetscScalar *array = m_raw.array;

        ScalarView values(array, nnz);
        IndexView row_ptr(ia, n + 1);
        IndexView colidx(ja, nnz);

        CRSMatrix<ScalarView, IndexView> crs(row_ptr, colidx, values, n);
        ilu_decompose(crs, modified);
    }

    void ILUDecompose<PetscMatrix, PETSC>::apply(const PetscMatrix &ilu, const PetscVector &b, PetscVector &x) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscVector L_inv_b(layout(b), 0.0);
        x.set(0.0);

        auto b_view = const_local_view_device(b);
        auto L_inv_b_view = local_view_device(L_inv_b);
        auto x_view = local_view_device(x);

        PetscSeqAIJRaw m_raw(ilu.raw_type());

        PetscInt n = m_raw.n;
        PetscInt nnz = m_raw.ia[n];

        const PetscInt *ia = m_raw.ia;
        const PetscInt *ja = m_raw.ja;
        PetscScalar *array = m_raw.array;

        ScalarView values(array, nnz);
        IndexView row_ptr(ia, n + 1);
        IndexView colidx(ja, nnz);

        CRSMatrix<ScalarView, IndexView> crs(row_ptr, colidx, values, n);

        auto b_array = b_view.array();
        auto L_array = L_inv_b_view.array();
        auto x_array = x_view.array();

        ilu_apply(crs, b_array, L_array, x_array);
    }

    void ILUDecompose<PetscMatrix, PETSC>::block_decompose(const PetscMatrix &, PetscMatrix &, const bool) {}

}  // namespace utopia
