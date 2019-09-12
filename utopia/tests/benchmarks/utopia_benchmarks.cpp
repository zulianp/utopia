#include "utopia_benchmarks.hpp"
#include "utopia_Benchmark_BLAS1.hpp"
#include "utopia_Benchmark_BLAS2.hpp"
#include "utopia_Benchmark_BLAS3.hpp"
#include "utopia_Benchmark_Algorithms.hpp"
#include "utopia_Benchmark_Access.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    void run_all_benchmarks(const std::string &backend_name)
    {
        if(mpi_world_rank() == 0) {
            std::cout << "\n> " << backend_name << std::endl;
        }

        int verbosity_level = 1;
        if(Utopia::instance().verbose()) {
            verbosity_level = 2;
        }

        BenchmarkBlas1<Matrix, Vector> blas1;
        blas1.set_verbosity_level(verbosity_level);
        blas1.run();

        BenchmarkBlas2<Matrix, Vector> blas2;
        blas2.set_verbosity_level(verbosity_level);
        blas2.run();

        BenchmarkBlas3<Matrix, Vector> blas3;
        blas3.set_verbosity_level(verbosity_level);
        blas3.run();

        BenchmarkAlgorithms<Matrix, Vector> algorithms;
        algorithms.set_verbosity_level(verbosity_level);
        algorithms.run();

        BenchmarkAccess<Matrix, Vector> access;
        access.set_verbosity_level(verbosity_level);
        access.run();
    }

    void run_benchmarks()
    {
        mpi_world_barrier();
        if(mpi_world_rank() == 0) {
            std::cout << "[Begin benchmark]" << std::endl;
        }

        //Parallel benchmarks
#ifdef WITH_PETSC
        run_all_benchmarks<DSMatrixd, DVectord>("petsc");
#ifdef PETSC_HAVE_CUDA
        run_all_benchmarks<CuSMatrixd, CuVectord>("petsc+cuda");
#endif //PETSC_HAVE_CUDA
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        run_all_benchmarks<TSMatrixd, TVectord>("trilinos");
#endif //WITH_TRILINOS

        //Serial benchmarks
        //FIXME does not compile for blas3 (missing mat*mat)
#ifdef WITH_BLAS
		if(mpi_world_size() == 1) {
			run_all_benchmarks<Matrixd, Vectord>("homemade");
		}
#endif

        mpi_world_barrier();
        if(mpi_world_rank() == 0) {
            std::cout << "[End benchmark]" << std::endl;
        }
    }
}