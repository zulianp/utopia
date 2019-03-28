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
        assemble_laplacian_1D(m);
        // disp(m);

        CuVectord x = local_values(ls.get(0), 1.);
        CuVectord y;
        y = m * x;

        CuVectord sol = local_zeros(local_size(y).get(0));

        //test mat-residual
        ConjugateGradient<CuSMatrixd, CuVectord, HOMEMADE>  cg;
        cg.max_it(n);
        cg.solve(m, y, sol);

        CuVectord r = y - m * sol;
        double err = norm2(r);
        utopia_test_assert(err < 1e-10);
    }

    template<class Matrix, class Vector>
    static void generic_linear_solver_test()
    {
        // SizeType n = 1000;
        SizeType n = 10000;
        Size ls{ n, n };
        Size gs{ ls.get(0) * mpi_world_size(), ls.get(1) * mpi_world_size() };

        Matrix m = local_sparse(ls.get(0), ls.get(1), 3);
        assemble_laplacian_1D(m);

        Vector x = local_values(ls.get(0), 0.);
        Vector b = local_values(ls.get(0), 1.);

        //test mat-residual
        KSPSolver<Matrix, Vector> ksp;
        // ksp.pc_type("sor");
        ksp.pc_type("jacobi");
        ksp.ksp_type("cg");
        ksp.max_it(n);
        ksp.solve(m, b, x);

        Vector r = b - m * x;
        double err = norm2(r);
        utopia_test_assert(err < 1e-10);
    }

    void petsc_cuda_linear_solver()
    {
       generic_linear_solver_test<CuSMatrixd, CuVectord>();
    }

    void petsc_no_cuda_linear_solver()
    {
        //normal matrices for comparing performance
        generic_linear_solver_test<DSMatrixd, DVectord>();
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

        multigrid.set_transfer_operators(std::move(interpolation_operators));
        multigrid.max_it(20);
        multigrid.atol(1e-10);
        multigrid.stol(1e-10);
        multigrid.rtol(1e-10);
  //      multigrid.verbose(true);

        CuVectord x = zeros(A.size().get(0));
        multigrid.update(make_ref(A));

        multigrid.apply(rhs, x);

        const double err = norm2(A * x - rhs);
        utopia_test_assert(err < 1e-6);
    }

#endif //PETSC_HAVE_CUDA


    void run_petsc_cuda_test() {
#ifdef PETSC_HAVE_CUDA
        UTOPIA_UNIT_TEST_BEGIN("PetscCudaTest");
        // UTOPIA_RUN_TEST(petsc_cuda_init);
        // UTOPIA_RUN_TEST(petsc_cuda_mg);
        UTOPIA_RUN_TEST(petsc_cuda_linear_solver);
        UTOPIA_RUN_TEST(petsc_no_cuda_linear_solver);
        UTOPIA_UNIT_TEST_END("PetscCudaTest");
#endif //PETSC_HAVE_CUDA
    }
}
