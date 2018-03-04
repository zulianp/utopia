
#include "utopia_LapackEigenSolver.hpp"
#include "utopia_blas_Backend.hpp"

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

	void spteqr_(char * compz, int * n, float * D, float * E, float * Z, int * ldz, float * work, int * info);
	void ssytrd_(char *uplo, int * n, float *A, int * lda, float * D, float * E, float * TAU, float * work, int * lwork, int * info);
	void sorgtr_(char *uplo, int * n, float *A, int * lda, float * tau, float * work, int *lwork, int *info);

}

namespace utopia {

	bool LapackEigenSolver::spd_geig(const Matrixd &A, const Matrixd &B, Vectord &evalues, Matrixd &evectors) const
	{

		const Size sA = size(A);
		const Size sB = size(B);

		Matrixd Acopy = A;
		Matrixd Bcopy = B;

		raw_type(evalues).resize(sA.get(0));
		raw_type(evectors).resize(sA.get(0), sB.get(0));

		int itype	= 1;
		char jobz 	= 'V';
		char range = 'A';
		char uplo 	= 'U';
		int n 		= sA.get(0);
		double * a = raw_type(Acopy).ptr();
		int lda	= n;
		double *b  = raw_type(Bcopy).ptr();
		int ldb   	= sB.get(0);
		double vl  = -1.0;
		double vu  = -1.0;
		int il 	= -1;
		int iu 	= -1;
		char cmach = 'S';
		double abstol =  2. * dlamch_(&cmach);
		int m 		= -1;
		double *w 	= &raw_type(evalues)[0];
		double *z 	= raw_type(evectors).ptr();
		int ldz 	= n;
		std::vector<double> work(8 * n, 0.0);
		int lwork = work.size();
		std::vector<int> iwork(5*n);
		std::vector<int> ifail(n);
		int info = -1;

		dsygvx_(&itype,
			&jobz,
			&range,
			&uplo,
			&n,
			a,
			&lda,
			b,
			&ldb,
			&vl,
			&vu,
			&il,
			&iu,
			&abstol,
			&m,
			w,
			z,
			&ldz,
			&work[0],
			&lwork,
			&iwork[0],
			&ifail[0],
			&info);

		if(info < 0) {
			std::cerr << "dsygvx_ " << (-info) << " argument is illegal" << std::endl;
		} else if(info > 0) {
			if(info > n) {
				std::cerr << "factorization could not be completed. Reason " << (info - n) << std::endl;
			} else {
				std::cerr << "some eigen vectors failed to converge" << std::endl;
				for(auto i : ifail) {
					std::cerr << i << ", ";
				}

				std::cerr << "\n";
			}
		}

		if(info == 0) return fix_sizes(m, evalues, evectors);
		return false;
	}

	bool LapackEigenSolver::spd_geig_small(const Matrixd &A, const Matrixd &B, const double upper_bound, Vectord &evalues, Matrixd &evectors) const
	{

		const Size sA = size(A);
		const Size sB = size(B);

		Matrixd Acopy = A;
		Matrixd Bcopy = B;

		raw_type(evalues).resize(sA.get(0));
		raw_type(evectors).resize(sA.get(0), sB.get(0));

		int itype	= 1;
		char jobz 	= 'V';
		char range	= 'V';
		char uplo 	= 'U';
		int n 		= sA.get(0);
		double * a 	= raw_type(Acopy).ptr();
		int lda		= n;
		double *b  	= raw_type(Bcopy).ptr();
		int ldb   	= sB.get(0);
		double vl  	= -1.0;
		double vu  	= upper_bound;
		int il 		= -1;
		int iu 		= -1;
		char cmach 	= 'S';
		double abstol =  2. * dlamch_(&cmach);
		int m 		= -1;
		double *w 	= &raw_type(evalues)[0];
		double *z 	=  raw_type(evectors).ptr();
		int ldz 	= n;
		std::vector<double> work(8 * n, 0.0);
		int lwork = work.size();
		std::vector<int> iwork(5*n);
		std::vector<int> ifail(n);
		int info = -1;

		dsygvx_(&itype,
			&jobz,
			&range,
			&uplo,
			&n,
			a,
			&lda,
			b,
			&ldb,
			&vl,
			&vu,
			&il,
			&iu,
			&abstol,
			&m,
			w,
			z,
			&ldz,
			&work[0],
			&lwork,
			&iwork[0],
			&ifail[0],
			&info);

		if(info < 0) {
			std::cerr << "dsygvx_: " << (-info) << " argument is illegal" << std::endl;
		} else if(info > 0) {
			if(info > n) {
				std::cerr << "dsygvx_: factorization could not be completed. Reason " << (info - n) << std::endl;
			} else {
				std::cerr << "dsygvx_: some eigen vectors failed to converge" << std::endl;
				for(auto i : ifail) {
					std::cerr << i << ", ";
				}

				std::cerr << "\n";
			}
		}

		if(info == 0) return fix_sizes(m, evalues, evectors);
		return false;
	}

	bool LapackEigenSolver::fix_sizes(const SizeType m, Vectord &evalues, Matrixd &evectors) const
	{
		const SizeType n = size(evalues).get(0);

		if(m == 0) {
			std::cerr << "wrong size: " <<  m << std::endl;
			return false;
		} else if(m != n) {
			evalues  = evalues.range(0, m);
			evectors = evectors.range(0, n, 0, m);
			return true;
		} 

		return true;
	}


	static bool aux_tridiagonalize(
		int n,
		const std::vector<double> &mat,
		std::vector<double> &diagonal,
		std::vector<double> &off_diagonal,
		std::vector<double> &reflectors
		)
	{
		char uplo = 'L';
		int lda = n;
		int lwork = -1;
		int info = -1;

		reflectors = mat;

		diagonal.resize(n);
		off_diagonal.resize(n-1);

		std::vector<double> tau(n -1);
		std::vector<double> work(1);

			//guess work size
		dsytrd_(
			&uplo,
			&n,
			&reflectors[0],
			&lda,
			&diagonal[0],
			&off_diagonal[0],
			&tau[0],
			&work[0],
			&lwork,
			&info);

		assert(info == 0);
		if(info != 0) return false;

		lwork = work[0];

		if(lwork > work.size()) {
			work.resize(lwork);
		}

		dsytrd_(
			&uplo,
			&n,
			&reflectors[0],
			&lda,
			&diagonal[0],
			&off_diagonal[0],
			&tau[0],
			&work[0],
			&lwork,
			&info);

		assert(info == 0);
		if(info != 0) return false;

		lwork = -1;
		dorgtr_(
			&uplo,
			&n,
			&reflectors[0],
			&lda, 
			&tau[0],
			&work[0], 
			&lwork,
			&info);

		assert(info == 0);

		lwork = work[0];

		if(lwork > work.size()) {
			work.resize(lwork);
		}

		dorgtr_(
			&uplo,
			&n,
			&reflectors[0],
			&lda, 
			&tau[0],
			&work[0], 
			&lwork,
			&info);

		assert(info == 0);
		return (info == 0);
	}

	static bool aux_eig(
		int n,
		const std::vector<double> &mat,
		std::vector<double> &eigen_values,
		std::vector<double> &eigen_vectors
		)
	{
		std::vector<double> work(n * 4);
		std::vector<double> off_diag(n - 1);

		char compz = 'V';
		int ldz = n;
		int info = 0;

		eigen_values.resize(n);
		eigen_vectors.resize(n*n);

		std::fill(eigen_vectors.begin(), eigen_vectors.end(), -6.);

		if(!aux_tridiagonalize(n, mat, eigen_values, off_diag, eigen_vectors)) 
		{
			return false;
		}

		dpteqr_(
			&compz,
			&n,
			&eigen_values[0],
			&off_diag[0],
			&eigen_vectors[0],
			&ldz,
			&work[0],
			&info);

		if(info != 0) {
			std::cerr << "Error in eig " << info << std::endl;
		} 

		for(int i = 0; i < n; ++i) {
			int offset_i = i * n;
			for(int j = i+1; j < n; ++j) {
				int offset_j = j * n;

				std::swap(eigen_vectors[offset_j + i],eigen_vectors[offset_i + j]);
			}
		}

		return (info == 0);
	}


	static bool aux_tridiagonalize(
		int n,
		const std::vector<float> &mat,
		std::vector<float> &diagonal,
		std::vector<float> &off_diagonal,
		std::vector<float> &reflectors
		)
	{
		char uplo = 'L';
		int lda = n;
		int lwork = -1;
		int info = -1;

		reflectors = mat;

		diagonal.resize(n);
		off_diagonal.resize(n-1);

		std::vector<float> tau(n -1);
		std::vector<float> work(1);

				//guess work size
		ssytrd_(
			&uplo,
			&n,
			&reflectors[0],
			&lda,
			&diagonal[0],
			&off_diagonal[0],
			&tau[0],
			&work[0],
			&lwork,
			&info);

		assert(info == 0);
		if(info != 0) return false;

		lwork = work[0];

		if(lwork > work.size()) {
			work.resize(lwork);
		}

		ssytrd_(
			&uplo,
			&n,
			&reflectors[0],
			&lda,
			&diagonal[0],
			&off_diagonal[0],
			&tau[0],
			&work[0],
			&lwork,
			&info);

		assert(info == 0);
		if(info != 0) return false;

		lwork = -1;
		sorgtr_(
			&uplo,
			&n,
			&reflectors[0],
			&lda, 
			&tau[0],
			&work[0], 
			&lwork,
			&info);

		assert(info == 0);
		if(info != 0) return false;

		lwork = work[0];

		if(lwork > work.size()) {
			work.resize(lwork);
		}

		sorgtr_(
			&uplo,
			&n,
			&reflectors[0],
			&lda, 
			&tau[0],
			&work[0], 
			&lwork,
			&info);

		assert(info == 0);
		return info == 0;
	}

	static bool aux_eig(
		int n,
		const std::vector<float> &mat,
		std::vector<float> &eigen_values,
		std::vector<float> &eigen_vectors
		)
	{
		std::vector<float> work(n * 4);
		std::vector<float> off_diag(n - 1);

		char compz = 'V';
		int ldz = n;
		int info = 0;

		eigen_values.resize(n);
		eigen_vectors.resize(n*n);

		std::fill(eigen_vectors.begin(), eigen_vectors.end(), -6.);

		aux_tridiagonalize(n, mat, eigen_values, off_diag, eigen_vectors);

		spteqr_(
			&compz,
			&n,
			&eigen_values[0],
			&off_diag[0],
			&eigen_vectors[0],
			&ldz,
			&work[0],
			&info);

		if(info != 0) {
			std::cerr << "Error in eig " << info << std::endl;
			return false;
		} 

		for(int i = 0; i < n; ++i) {
			int offset_i = i * n;
			for(int j = i+1; j < n; ++j) {
				int offset_j = j * n;

				std::swap(eigen_vectors[offset_j + i], eigen_vectors[offset_i + j]);
			}
		}

		return true;
	}

	bool LapackEigenSolver::spd_eig(const Matrixd &A, Vectord &evalues, Matrixd &evector) const
	{
		auto n = size(A).get(0);

		evector = zeros(size(A));
		evalues = zeros(n);

		return aux_eig(n, A.implementation().entries(), evalues.implementation(), evector.implementation().entries());
	}
}
