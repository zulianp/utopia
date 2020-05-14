#include "utopia_lapack_Algorithms.hpp"


extern "C" {
    int dsygvx_(int *itype,
        char *jobz,
        char *range,
        char * uplo,
        int *n,
        double *a,
        int *lda,
        double *b,
        int *ldb,
        double *vl,
        double *vu,
        int *il,
        int *iu,
        double *abstol,
        int *m,
        double *w,
        double *z__,
        int *ldz,
        double *work,
        int *lwork,
        int *iwork,
        int *ifail,
        int *info);

    double dlamch_(char *cmach);


    void dpteqr_(char * compz, int * n, double * D, double * E, double * Z, int * ldz, double * work, int * info);
    void dsytrd_(char *uplo, int * n, double *A, int * lda, double * D, double * E, double * TAU, double * work, int * lwork, int * info);
    void dorgtr_(char *uplo, int * n, double *A, int * lda, double * tau, double * work, int *lwork, int *info);

    void dsyev_(
        char * jobz,
        char * uplo,
        int * n,
        double * a,
        int * lda,
        double * w,
        double * work,
        int * lwork,
        int * info
    );

    // void spteqr_(char * compz, int * n, float * D, float * E, float * Z, int * ldz, float * work, int * info);
    // void ssytrd_(char *uplo, int * n, float *A, int * lda, float * D, float * E, float * TAU, float * work, int * lwork, int * info);
    // void sorgtr_(char *uplo, int * n, float *A, int * lda, float * tau, float * work, int *lwork, int *info);

}

namespace utopia {

    template<>
    int LapackAlgorithms<double>::sygvx(
        int itype,
        char jobz,
        char range,
        char  uplo,
        int n,
        double * a,
        int lda,
        double * b,
        int ldb,
        double * vl,
        double * vu,
        int il,
        int iu,
        double abstol,
        int m,
        double * w,
        double * z,
        int ldz,
        double * work,
        int lwork,
        int iwork,
        int ifail,
        int *info)
    {
        return dsygvx_(
            &itype,
            &jobz,
            &range,
            &uplo,
            &n,
            a,
            &lda,
            b,
            &ldb,
            vl,
            vu,
            &il,
            &iu,
            &abstol,
            &m,
            w,
            z,
            &ldz,
            work,
            &lwork,
            &iwork,
            &ifail,
            info
        );
    }

    template<>
    double LapackAlgorithms<double>::lamch(char cmach)
    {
        return dlamch_(&cmach);
    }

    template<>
    void LapackAlgorithms<double>::pteqr(char compz, int  n, double * D, double * E, double * Z, int  ldz, double * work, int *info)
    {
        dpteqr_(
            &compz,
            &n,
             D,
             E,
             Z,
            &ldz,
             work,
            info
        );
    }

    template<>
    void LapackAlgorithms<double>::sytrd(char uplo, int  n, double * A, int lda, double * D, double * E, double * TAU, double * work, int lwork, int *info)
    {
        dsytrd_(
            &uplo,
            &n,
            A,
            &lda,
            D,
            E,
            TAU,
            work,
            &lwork,
            info
        );
    }

    template<>
    void LapackAlgorithms<double>::orgtr(char uplo, int  n, double * A, int lda, double * tau, double * work, int lwork, int *info)
    {
        dorgtr_(
            &uplo,
            &n,
            A,
            &lda,
            tau,
            work,
            &lwork,
            info
        );
    }

    template<>
    void LapackAlgorithms<double>::syev(
        char jobz,
        char uplo,
        int n,
        double * a,
        int lda,
        double * w,
        double * work,
        int lwork,
        int * info
    )
    {
        dsyev_(
            &jobz,
            &uplo,
            &n,
            a,
            &lda,
            w,
            work,
            &lwork,
            info
        );

    }

}  // namespace utopia
