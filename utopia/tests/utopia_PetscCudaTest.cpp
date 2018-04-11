#include "utopia.hpp"
#include "utopia_PetscCudaTest.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

namespace utopia {

#ifdef WITH_PETSC_CUDA
    void petsc_cuda_init()
    {
    	//DSMatrixd m = local_sparse(10, 10, 1, sub_comm, str("my_mat"), device::gpu);

    	Size ls{10, 10};
    	Size gs{ls.get(0) * mpi_world_size(), ls.get(1) * mpi_world_size() };

    	DSMatrixd m;
    	m.implementation().mat_aij_cusparse_init(
    		PETSC_COMM_WORLD, ls.get(0), ls.get(1), gs.get(0), gs.get(1), 3, 3
    	);

    	assemble_laplacian_1D(size(m).get(0), m);

    	disp(m);

    	DVectord x;
    	x.implementation().values(
    		PETSC_COMM_WORLD, VECMPICUDA, ls.get(1), gs.get(1), 1.
    	);


    	DVectord y;
    	y.implementation().values(
    		PETSC_COMM_WORLD, VECMPICUDA, ls.get(0), gs.get(0), 2.
    	);

    	y = m * x;

    	disp(y);
    }

#endif //WITH_PETSC_CUDA;

    void run_petsc_cuda_test() {
#ifdef WITH_PETSC_CUDA
        UTOPIA_UNIT_TEST_BEGIN("PetscCudaTest");
        UTOPIA_RUN_TEST(petsc_cuda_init);
        UTOPIA_UNIT_TEST_END("PetscCudaTest");
#endif // WITH_PETSC_CUDA
    }
}
