#include "utopia.hpp"
#include "utopia_PetscCudaTest.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

namespace utopia {

#ifdef PETSC_HAVE_CUDA

    void petsc_cuda_init()
    {
        SizeType n = 1000;
        Size ls{ n, n };
        Size gs{ ls.get(0) * mpi_world_size(), ls.get(1) * mpi_world_size() };

        CuSMatrixd m = local_sparse(ls.get(0), ls.get(1), 3);
        assemble_laplacian_1D(size(m).get(0), m);
        // disp(m);

        CuVectord x = local_values(ls.get(0), 1.);
        CuVectord y;
        y = m * x;

        // disp(y);

        CuVectord sol = local_zeros(local_size(y).get(0));

        //test mat-residual
        ConjugateGradient<CuSMatrixd, CuVectord, HOMEMADE>  cg;
        cg.verbose(true);
        cg.max_it(n);
        cg.solve(m, y, sol);

        CuVectord r = y - m * sol; 
        double err = norm2(r);
        assert(err < 1e-10);
    }

    void petsc_cuda_mg()
    {
        CuVectord rhs;
        CuSMatrixd A, I_1, I_2, I_3;

        const std::string data_path = Utopia::instance().get("data_path");

        read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
        read(data_path + "/laplace/matrices_for_petsc/f_A", A);
        read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
        read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

        std::vector<std::shared_ptr<CuSMatrixd>> interpolation_operators;
        interpolation_operators.push_back(make_ref(I_2));
        interpolation_operators.push_back(make_ref(I_3));

        auto smoother      = std::make_shared<ConjugateGradient<CuSMatrixd, CuVectord, HOMEMADE>>();
        auto linear_solver = std::make_shared<ConjugateGradient<CuSMatrixd, CuVectord, HOMEMADE>>();

        Multigrid<CuSMatrixd, CuVectord> multigrid(smoother, linear_solver);

        // smoother->verbose(true);
        // linear_solver->verbose(true);

        multigrid.init_transfer_from_fine_to_coarse(std::move(interpolation_operators));
        multigrid.max_it(20);
        multigrid.atol(1e-10);
        multigrid.stol(1e-10);
        multigrid.rtol(1e-10);
        multigrid.verbose(true);

        CuVectord x = zeros(A.size().get(0));
        multigrid.update(make_ref(A));

        multigrid.apply(rhs, x);

        const double err = norm2(A * x - rhs);
        assert(err < 1e-6);
    }

#endif //PETSC_HAVE_CUDA


    void run_petsc_cuda_test() {
#ifdef PETSC_HAVE_CUDA
        UTOPIA_UNIT_TEST_BEGIN("PetscCudaTest");
        UTOPIA_RUN_TEST(petsc_cuda_init);
        UTOPIA_RUN_TEST(petsc_cuda_mg);
        UTOPIA_UNIT_TEST_END("PetscCudaTest");
#endif //PETSC_HAVE_CUDA
    }
}
