#include "utopia_petsc_Eval_Multiply.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Logger.hpp"
#include "utopia_petsc_Matrix.hpp"

#include "utopia_Tracer.hpp"

namespace utopia {

    inline static bool check_error(const SizeType err) { return PetscErrorHandler::Check(err); }

    void PetscEvalTripleMatrixProduct::ptap(PetscMatrix &result,
                                            const PetscMatrix &A,
                                            const PetscMatrix &P,
                                            MatReuse reuse) {
        UTOPIA_TRACE_SCOPE("PetscEvalTripleMatrixProduct::ptap");

        if (result.empty()) {
            // Safegard
            reuse = MAT_INITIAL_MATRIX;
        }

        if (reuse == MAT_REUSE_MATRIX) {
            check_error(MatPtAP(raw_type(A), raw_type(P), MAT_REUSE_MATRIX, 1., &raw_type(result)));
            return;
        }

        if (!A.is_cuda() && !A.is_block() && (result.is_alias(A) || result.is_alias(P))) {
            Mat temp = result.raw_type();

            check_error(MatPtAP(A.raw_type(), P.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()));

            MatDestroy(&temp);
        } else {
            MatDestroy(&result.raw_type());

            if (A.is_cuda()) {
                m_utopia_status_once(
                    "MatPtAP does not work properly with the cusparse backend. Workaround implemented.");
                Mat temp;
                check_error(MatPtAP(A.raw_type(), P.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp));
                check_error(MatConvert(temp, A.type(), MAT_INITIAL_MATRIX, &result.raw_type()));
                MatDestroy(&temp);
            } else if (A.is_block()) {
                m_utopia_status_once("MatPtAP does not work with the matbaij type. Workaround with copy implemented.");

                Mat temp;
                check_error(MatConvert(A.raw_type(), P.type(), MAT_INITIAL_MATRIX, &temp));

                check_error(MatPtAP(temp, P.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp));
                check_error(MatConvert(temp, A.type(), MAT_INITIAL_MATRIX, &result.raw_type()));

                MatDestroy(&temp);

            } else {
                check_error(MatPtAP(A.raw_type(), P.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()));
            }
        }

        assert(result.same_type(A));
    }

    void PetscEvalTripleMatrixProduct::rart(PetscMatrix &result,
                                            const PetscMatrix &A,
                                            const PetscMatrix &R,
                                            MatReuse reuse) {
        if (result.empty()) {
            // Safegard
            reuse = MAT_INITIAL_MATRIX;
        }

        if (result.is_alias(A) || result.is_alias(R)) {
            Mat temp = result.raw_type();

            check_error(MatRARt(A.raw_type(), R.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()));

            MatDestroy(&temp);
        } else {
            MatDestroy(&result.raw_type());
            check_error(MatRARt(A.raw_type(), R.raw_type(), reuse, PETSC_DEFAULT, &result.raw_type()));
        }
    }

    void PetscEvalTripleMatrixProduct::abc(PetscMatrix &result,
                                           const PetscMatrix &A,
                                           const PetscMatrix &B,
                                           const PetscMatrix &C,
                                           MatReuse reuse) {
        if (result.empty()) {
            // Safegard
            reuse = MAT_INITIAL_MATRIX;
        }

        if (result.is_alias(A) || result.is_alias(B) || result.is_alias(C)) {
            Mat temp = result.raw_type();

            check_error(MatMatMatMult(
                A.raw_type(), B.raw_type(), C.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()));

            MatDestroy(&temp);

        } else {
            result.destroy();
            check_error(
                MatMatMatMult(A.raw_type(), B.raw_type(), C.raw_type(), reuse, PETSC_DEFAULT, &result.raw_type()));
        }

        assert(result.same_type(A));
    }

    void ptap_reuse_matrix(const PetscMatrix &A, const PetscMatrix &P, PetscMatrix &result) {
        UTOPIA_TRACE_SCOPE("ptap_reuse_matrix");
        PetscEvalTripleMatrixProduct::ptap(result, A, P, MAT_REUSE_MATRIX);
    }

}  // namespace utopia
