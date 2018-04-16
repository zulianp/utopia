#include "utopia_benchmarks.hpp"
#include "utopia_Benchmark_BLAS1.hpp"

namespace utopia {
	void run_benchmarks()
	{
		mpi_world_barrier();
		if(mpi_world_rank() == 0) {
			std::cout << "[Begin benchmark]" << std::endl;
		}

		//Parallel benchmarks
#ifdef WITH_PETSC
		{
			if(mpi_world_rank() == 0) {
				std::cout << "> petsc" <<std::endl; 
			}

			BenchmarkBlas1<DSMatrixd, DVectord> petsc_b1;
			petsc_b1.run();
		}

#ifdef PETSC_HAVE_CUDA
		{
			if(mpi_world_rank() == 0) {
				std::cout << "> petsc+cuda" <<std::endl; 
			}

			BenchmarkBlas1<CuSMatrixd, CuVectord> petsc_b1;
			petsc_b1.run();
		}

#endif //PETSC_HAVE_CUDA

#endif //WITH_PETSC

#ifdef WITH_TRILINOS
		{
			if(mpi_world_rank() == 0) {
				std::cout << "> trilinos" <<std::endl; 
			}
			
			BenchmarkBlas1<TSMatrixd, TVectord> trilinos_b1;
			trilinos_b1.run();
		}

#endif //WITH_TRILINOS

		//Serial benchmarks
		if(mpi_world_size() == 1) {
			// BenchmarkBlas1<Matrixd, Vectord> blas_b1;
			// blas_b1.run();
		}

		mpi_world_barrier();
		if(mpi_world_rank() == 0) {
			std::cout << "[End benchmark]" << std::endl;
		}
	}
}