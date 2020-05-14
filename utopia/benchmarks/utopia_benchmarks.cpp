#include "utopia_benchmarks.hpp"
#include "utopia_Benchmark_Access.hpp"
#include "utopia_Benchmark_Algorithms.hpp"
#include "utopia_Benchmark_BLAS1.hpp"
#include "utopia_Benchmark_BLAS2.hpp"
#include "utopia_Benchmark_BLAS3.hpp"
#include "utopia_Benchmark_Transform.hpp"

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

        BenchmarkAccess<Matrix, Vector> access;
        access.set_verbosity_level(verbosity_level);
        access.run();

        BenchmarkTransform<Matrix, Vector> transform;
        transform.set_verbosity_level(verbosity_level);
        transform.run();

        BenchmarkAlgorithms<Matrix, Vector> algorithms;
        algorithms.set_verbosity_level(verbosity_level);
        algorithms.run();
    }

    void run_benchmarks()
    {
        mpi_world_barrier();
        if(mpi_world_rank() == 0) {
            std::cout << "[Begin benchmark]" << std::endl;
        }

        //Parallel benchmarks
#ifdef WITH_PETSC
        run_all_benchmarks<PetscMatrix, PetscVector>("petsc");
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        run_all_benchmarks<TpetraMatrixd, TpetraVectord>("trilinos");
#endif //WITH_TRILINOS

        //Serial benchmarks
#ifdef WITH_BLAS
		if(mpi_world_size() == 1) {
			run_all_benchmarks<BlasMatrixd, BlasVectord>("homemade");
		}
#endif

        mpi_world_barrier();
        if(mpi_world_rank() == 0) {
            std::cout << "[End benchmark]" << std::endl;
        }
    }
}  // namespace utopia
