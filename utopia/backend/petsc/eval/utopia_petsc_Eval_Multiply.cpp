#include "utopia_petsc_Eval_Multiply.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Instance.hpp"

namespace utopia {

    inline static bool check_error(const SizeType err) {
        return PetscErrorHandler::Check(err);
    }
    
    void PetscEvalTripleMatrixProduct::ptap(PetscMatrix & result, const PetscMatrix & A, const PetscMatrix &P)
    {
        if(result.raw_type() != A.raw_type() &&
           result.raw_type() != P.raw_type()) {
            MatDestroy(&result.raw_type());
        } else {
            std::cerr << "[Error] not handled case in ptap" << std::endl;
            assert(false);
        }

        assert(A.same_type(P));

        if(A.is_cuda()) {
            m_utopia_status_once("MatPtAP does not work properly with the cusparse backend. Workaround implemented.");
            Mat temp;
            check_error( MatPtAP(A.raw_type(), P.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp) );
            check_error( MatConvert(temp, A.type(), MAT_INITIAL_MATRIX, &result.raw_type()) );
            MatDestroy(&temp);
        } else {
            check_error( MatPtAP(A.raw_type(), P.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()) );
        }

        assert(result.same_type(A));
    }

    void PetscEvalTripleMatrixProduct::rart(PetscMatrix &result, const PetscMatrix &A, const PetscMatrix &R)
    {
        if(result.raw_type() != A.raw_type() &&
           result.raw_type() != R.raw_type()) {
            MatDestroy(&result.raw_type());
        } else {
            std::cerr << "[Error] not handled case in rart" << std::endl;
            assert(false);
        }

        check_error( MatRARt(A.raw_type(), R.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()) );
    }

    void PetscEvalTripleMatrixProduct::abc(PetscMatrix &result, const PetscMatrix &A, const PetscMatrix &B, const PetscMatrix &C)
    {
        if(result.raw_type() != A.raw_type() && result.raw_type() != B.raw_type() && result.raw_type() != C.raw_type()) {
            result.destroy();
        }

        check_error( MatMatMatMult(A.raw_type(), B.raw_type(), C.raw_type(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result.raw_type()) );

        assert(result.same_type(A));
    }

}