#ifndef UTOPIA_BLAS_ADAPTER_HPP
#define UTOPIA_BLAS_ADAPTER_HPP

namespace utopia {
    template <typename T>
    class BLASAlgorithms {
    public:
        static void gemm(const char transa,
                         const char transb,
                         const int m,
                         const int n,
                         const int k,
                         const T alpha,
                         const T *a,
                         const int lda,
                         const T *b,
                         const int ldb,
                         const T beta,
                         T *c,
                         const int ldc);

        static void gemv(const char trans,
                         const int m,
                         const int n,
                         const T alpha,
                         const T *a,
                         const int lda,
                         const T *x,
                         const int incx,
                         const T beta,
                         T *c,
                         const int incy);

        static void axpy(const int n, const T alpha, const T *x, const int incx, T *y, const int incy);

        static T nrm2(const int n, const T *x, const int incx);
        static T asum(const int n, const T *x, const int incx);
        static int amax(const int n, const T *x, const int incx);

        static T ddot(const int n, const T *sx, const int incx, const T *sy, const int incy);

        static void scal(const int n, const T a, T *sx, const int incx);
        static void copy(const int n, const double *x, const int incx, double *y, const int incy);
    };

}  // namespace utopia

#endif  // UTOPIA_BLAS_ADAPTER_HPP
