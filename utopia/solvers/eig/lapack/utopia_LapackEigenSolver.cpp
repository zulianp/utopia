
#include "utopia_LapackEigenSolver.hpp"

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
}