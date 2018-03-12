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
		BenchmarkBlas1<DSMatrixd, DVectord> petsc_b1;
		petsc_b1.run();


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