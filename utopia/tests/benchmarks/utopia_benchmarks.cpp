#include "utopia_benchmarks.hpp"
#include "utopia_Benchmark_BLAS1.hpp"

namespace utopia {

	template<class Matrix, class Vector>
	void run_all_benchmarks(const std::string &backend_name)
	{
		if(mpi_world_rank() == 0) {
			std::cout << "> " << backend_name <<std::endl; 
		}

		BenchmarkBlas1<Matrix, Vector> blas1;
		blas1.run();
	}

	void run_benchmarks()
	{
		mpi_world_barrier();
		if(mpi_world_rank() == 0) {
			std::cout << "[Begin benchmark]" << std::endl;
		}

		//Parallel benchmarks
#ifdef WITH_PETSC
		run_all_benchmarks<DSMatrixd, Vectord>("petsc");
#ifdef PETSC_HAVE_CUDA
		run_all_benchmarks<CuSMatrixd, CuVectord>("petsc+cuda");
#endif //PETSC_HAVE_CUDA
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
		run_all_benchmarks<TSMatrixd, TVectord>("trilinos");
#endif //WITH_TRILINOS

		//Serial benchmarks
#ifdef WITH_BLAS
		if(mpi_world_size() == 1) {
			run_all_benchmarks<CRSMatrixd, Vectord>("homemade");
		}
#endif	

		mpi_world_barrier();
		if(mpi_world_rank() == 0) {
			std::cout << "[End benchmark]" << std::endl;
		}
	}
}