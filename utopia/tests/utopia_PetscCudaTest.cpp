#include "utopia.hpp"
#include "utopia_PetscCudaTest.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

namespace utopia {

//#ifdef WITH_PETSC_CUDA
    void petsc_cuda_init()
    {
    	Size ls{ 10, 10 };
    	Size gs{ ls.get(0) * mpi_world_size(), ls.get(1) * mpi_world_size() };

    	CuSMatrixd m = local_sparse(ls.get(0), ls.get(1), 3);
    	assemble_laplacian_1D(size(m).get(0), m);
    	disp(m);

    	CuVectord x = local_values(ls.get(0), 1.);
    	CuVectord y;
    	y = m * x;

    	disp(y);

        CuVectord sol = local_zeros(local_size(y).get(0));

        //test mat-redisual
        CuVectord r = y - m * sol; 

        // ConjugateGradient<CuSMatrixd, CuVectord, HOMEMADE>  cg;
        // cg.verbose(true);
        // cg.solve(m, y, sol);

        disp(r);
    }

//#endif //WITH_PETSC_CUDA;

    void run_petsc_cuda_test() {
//#ifdef WITH_PETSC_CUDA
        UTOPIA_UNIT_TEST_BEGIN("PetscCudaTest");
        UTOPIA_RUN_TEST(petsc_cuda_init);
        UTOPIA_UNIT_TEST_END("PetscCudaTest");
//#endif // WITH_PETSC_CUDA
    }
}
