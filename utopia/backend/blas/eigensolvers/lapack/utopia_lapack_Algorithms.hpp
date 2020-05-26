#ifndef UTOPIA_LAPACK_ALGORITHMS_HPP
#define UTOPIA_LAPACK_ALGORITHMS_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Operators.hpp"

namespace utopia {
    template <typename T>
    class LapackAlgorithms {
    public:
        static int sygvx(int itype,
                         char jobz,
                         char range,
                         char uplo,
                         int n,
                         T *a,
                         int lda,
                         T *b,
                         int ldb,
                         T *vl,
                         T *vu,
                         int il,
                         int iu,
                         T abstol,
                         int m,
                         T *w,
                         T *z,
                         int ldz,
                         T *work,
                         int lwork,
                         int iwork,
                         int ifail,
                         int *info);

        static T lamch(char cmach);
        static void pteqr(char compz, int n, T *D, T *E, T *Z, int ldz, T *work, int *info);
        static void sytrd(char uplo, int n, T *A, int lda, T *D, T *E, T *TAU, T *work, int lwork, int *info);
        static void orgtr(char uplo, int n, T *A, int lda, T *tau, T *work, int lwork, int *info);

        static void
        syev(char jobz, char uplo, int n, double *a, int lda, double *w, double *work, int lwork, int *info);
    };

    // template<typename T, int Dim>
    // class LapackEigenDecomposition {
    // public:

    //     bool apply(const T *A, const T*e, const T*v)
    //     {
    //         T work[Dim * 4];
    //         T off_diag[Dim - 1];
    //         T tau[Dim - 1];

    //         char uplo = 'L';
    //         int lda = Dim;
    //         int lwork = -1;
    //         int info = -1;
    //         int n = Dim;

    //         device::copy(A, A + Dim * Dim, v);

    //         syev(
    //             'V',
    //             char uplo,
    //             int n,
    //             double * a,
    //             int lda,
    //             double * w,
    //             double * work,
    //             int lwork,
    //             int * info
    //         );

    //     }

    // };

}  // namespace utopia

#endif  // UTOPIA_LAPACK_ALGORITHMS_HPP
