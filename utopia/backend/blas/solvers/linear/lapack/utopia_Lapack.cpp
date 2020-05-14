//
// Created by Alessandro Rigazzi on 22/05/15.
//

#include "utopia_Lapack.hpp"

extern "C"
{
    void dgesv_(const int* n, const int* nrhs, const double* a, const int* lda,
     const int* ipiv, const double* b, const int* ldb, const int* info );
}

namespace utopia {
    namespace internals {
        bool lapack_dgesv_solve(const BlasMatrixd&A, const BlasVectord&b, BlasVectord&x)
        {
            const int n = A.rows();
            const int nrhs = 1;
            const int lda = n;
            const int ldb = n;
            int ipiv[n];
            int info;

            x = b;

            BlasMatrixd Atmp = A;

            dgesv_(&n, &nrhs, Atmp.ptr(), &lda, ipiv, x.ptr(), &ldb, &info);
            return true;
        }
    }  // namespace internals
}  // namespace utopia
