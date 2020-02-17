#include "utopia_blas_Algorithms.hpp"

//Try: https://github.com/flame/blis#documentation

extern "C" {
    // Fortran interface (taken from www.netlib.org/clapack)...
    void dgemm_(const char *transa, const char *transb, const int *m,
                const int *n, const int *k,
                const double *alpha,
                const double *a, const int *lda, const double *b, const int *ldb,
                const double *beta,
                double *c, const int *ldc);

    void dgemv_(const char *trans, const int *m, const int *n,
                const double *alpha, const double *a, const int *lda, const double *x,
                const int *incx, const double *beta, double *c, const int *incy);


    void daxpy_(const int *n,
                const double *alpha,
                const double *x,
                const int *incx,
                double *y,
                const int *incy);

    double dnrm2_(const int *n, const double *x, const int *incx);
    double dasum_(const int *n, const double *x, const int *incx);
    int idamax_(const int *n, const double *x, const int *incx);

    double ddot_(const int *n, const double *sx, const int *incx, const double *sy, const int *incy);

    void dscal_(const int *n,
                const double *sa,
                double *sx,
                const int *incx);

    void dcopy_	(
    	const int *n,
    	const double * x,
    	const int *incx,
    	double * y,
    	const int *incy
    );


#ifdef WITH_MKL

    ///sparse-blas
    void mkl_dcsrgemv_(const char * transa,
                       const int * m,
                       const double * a,
                       const int * ia,
                       const int * ja,
                       const double * x,
                       double * y);

#endif //WITH_MKL

}

namespace utopia {

	template<>
	void BLASAlgorithms<double>::gemm(
		const char transa, const char transb, const int m,
		const int n, const int k,
		const double alpha,
		const double *a, const int lda, const double *b, const int ldb,
		const double beta,
		double *c, const int ldc)
	{
		dgemm_(
			&transa,
			&transb,
			&m,
			&n,
			&k,
			&alpha,
			a,
			&lda,
			b,
			&ldb,
			&beta,
		    c,
		    &ldc
		 );
	}

	template<>
	void BLASAlgorithms<double>::gemv(
		const char trans, const int m, const int n,
		const double alpha, const double *a, const int lda, const double *x,
		const int incx, const double beta, double *c, const int incy)
	{
		dgemv_(
			&trans,
			&m,
			&n,
		    &alpha,
		    a,
		    &lda,
		    x,
		    &incx,
		    &beta,
		    c,
		    &incy
		);
	}

	template<>
	void BLASAlgorithms<double>::axpy(
		const int n,
		const double alpha,
		const double *x,
		const int incx,
		double *y,
		const int incy)
	{
		daxpy_(
			&n,
			&alpha,
			x,
			&incx,
			y,
			&incy
		);
	}

	template<>
	double BLASAlgorithms<double>::nrm2(const int n, const double *x, const int incx)
	{
		return dnrm2_(
			&n,
			x,
			&incx
		);
	}

	template<>
	double BLASAlgorithms<double>::asum(const int n, const double *x, const int incx)
	{
		return dasum_(&n, x, &incx);
	}

	template<>
	int BLASAlgorithms<double>::amax(const int n, const double *x, const int incx)
	{
		return idamax_(&n, x, &incx);
	}

	template<>
	double BLASAlgorithms<double>::ddot(const int n, const double *sx, const int incx, const double *sy, const int incy)
	{
		return ddot_(&n, sx, &incx, sy, &incy);
	}

	template<>
	void BLASAlgorithms<double>::scal(const int n, const double a, double *sx, const int incx)
	{
		dscal_(&n, &a, sx, &incx);
	}

	template<>
	void BLASAlgorithms<double>::copy(const int n, const double *x, const int incx, double *y, const int incy)
	{
		dcopy_	(
			&n,
			x,
			&incx,
			y,
			&incy
		);
	}

	//explicit instantiation for double
	template class BLASAlgorithms<double>;
}
