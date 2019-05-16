#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_MultiLevelTestProblem.hpp"

namespace utopia {

#ifdef WITH_PETSC
    class PetscLinearSolverTest {
    public:

        void run()
        {
            UTOPIA_RUN_TEST(petsc_cg);
            UTOPIA_RUN_TEST(petsc_bicgstab);
            UTOPIA_RUN_TEST(petsc_gmres);
//FIXME make it work also without mumps
#ifdef PETSC_HAVE_MUMPS
            UTOPIA_RUN_TEST(petsc_mg);
            UTOPIA_RUN_TEST(petsc_cg_mg);
            UTOPIA_RUN_TEST(petsc_mg_1D);
            UTOPIA_RUN_TEST(petsc_block_mg_exp);
            UTOPIA_RUN_TEST(petsc_block_mg);
            UTOPIA_RUN_TEST(petsc_mg_exp);
            UTOPIA_RUN_TEST(petsc_superlu_mg);
            UTOPIA_RUN_TEST(petsc_mg_jacobi);
            UTOPIA_RUN_TEST(petsc_factorization);
            UTOPIA_RUN_TEST(petsc_st_cg_mg);

#endif //PETSC_HAVE_MUMPS
        }

        void petsc_cg()
        {
            MultiLevelTestProblem<DSMatrixd, DVectord> ml_problem(100, 2);
            DVectord x = zeros(size(*ml_problem.rhs));
            (*ml_problem.rhs) *= 0.0001;

            ConjugateGradient<DSMatrixd, DVectord, HOMEMADE> cg;
            cg.rtol(1e-6);
            cg.atol(1e-6);
            cg.max_it(500);
            // cg.verbose(true);

            x = *ml_problem.rhs;
            cg.update(ml_problem.matrix);
            cg.apply(*ml_problem.rhs, x);

            utopia_test_assert(approxeq(*ml_problem.rhs, *ml_problem.matrix * x, 1e-5));
        }

        void petsc_mg_1D()
        {
            if(mpi_world_size() > 1) return;

            const static bool verbose = false;

            MultiLevelTestProblem<DSMatrixd, DVectord> ml_problem(4, 2);
            // ml_problem.write_matlab("./");

            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            smoother->verbose(verbose);

            Multigrid<DSMatrixd, DVectord> multigrid(
                smoother,
                std::make_shared<Factorization<DSMatrixd, DVectord>>()
                // std::make_shared<ConjugateGradient<DSMatrixd, DVectord, HOMEMADE>>(),
                // std::make_shared<ConjugateGradient<DSMatrixd, DVectord, HOMEMADE>>()
            );

            multigrid.set_transfer_operators(ml_problem.interpolators);
            multigrid.max_it(4);
            multigrid.atol(1e-12);
            multigrid.stol(1e-10);
            multigrid.rtol(1e-10);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(verbose);

            DVectord x = zeros(size(*ml_problem.rhs));
            multigrid.update(ml_problem.matrix);

            if(verbose) {
                multigrid.describe();
            }

            multigrid.apply(*ml_problem.rhs, x);
            utopia_test_assert(approxeq(*ml_problem.rhs, *ml_problem.matrix * x, 1e-7));
        }

        void petsc_mg_exp()
        {

            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
            // auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
            // linear_solver->verbose(true);
            Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> multigrid;//(smoother, linear_solver);
            // multigrid.set_default_pc_type(PCILU);
            // multigrid.set_default_ksp_type(KSPFGMRES);

            // Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> multigrid;
            multigrid.set_transfer_operators(std::move(interpolation_operators));
            multigrid.max_it(20);
            multigrid.atol(1e-15);
            multigrid.stol(1e-15);
            multigrid.rtol(1e-15);

            DVectord x = zeros(A.size().get(0));
            // multigrid.verbose(true);
            multigrid.solve(A, rhs, x);

            const double err = norm2(A*x - rhs);
            utopia_test_assert(err < 1e-6);
        }

        void petsc_mg()
        {
            // reading data from outside
            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            // read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;

            // from coarse to fine
            // interpolation_operators.push_back(std::move(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            //  init
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(Solver::mumps(), Solver::lu_decomposition());
#endif //PETSC_HAVE_MUMPS

            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();

            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            // multigrid.set_use_line_search(true);
            // multigrid.verbose(true);


            multigrid.set_transfer_operators(std::move(interpolation_operators));
            // multigrid.set_fix_semidefinite_operators(true);
            multigrid.update(make_ref(A));

            DVectord x_0 = zeros(A.size().get(0));

            // multigrid.verbose(true);
            multigrid.apply(rhs, x_0);

            x_0 = zeros(A.size().get(0));
            multigrid.cycle_type(FULL_CYCLE);
            multigrid.apply(rhs, x_0);

            x_0 = zeros(A.size().get(0));
            multigrid.cycle_type(FULL_CYCLE);
            multigrid.v_cycle_repetition(2);

            multigrid.apply(rhs, x_0);

            multigrid.max_it(1);
            multigrid.cycle_type(MULTIPLICATIVE_CYCLE);
            auto gmres = std::make_shared<GMRES<DSMatrixd, DVectord>>();
            gmres->set_preconditioner(make_ref(multigrid));
            x_0.set(0.);
            gmres->verbose(false);
            gmres->atol(1e-16);
            gmres->rtol(1e-16);
            gmres->solve(A, rhs, x_0);

            double diff = norm2(rhs - A *x_0);
            if(diff > 1e-6) {
                utopia_error("petsc_mg: gmres preconditioned with mg does not do what it is supposed to");
            }
            // utopia_test_assert( approxeq(diff, 0., 1e-6) );
        }

        void petsc_bicgstab()
        {

            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);

            BiCGStab<DMatrixd, DVectord> bicgs;
            bicgs.solve(mat, rhs, sol);

            DVectord expected = zeros(_n);
            utopia_test_assert(approxeq(expected, sol));
        }

        void petsc_gmres()
        {
            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);

            GMRES<DMatrixd, DVectord> gmres;
            gmres.pc_type(PCBJACOBI);

            gmres.number_of_subdomains(mpi_world_size());
            gmres.update(std::make_shared<DMatrixd>(mat));
            gmres.sub_ksp_pc_type(KSPPREONLY, PCILU);

            gmres.verbose(false);
            gmres.solve(mat, rhs, sol);

            DVectord expected = zeros(_n);
            utopia_test_assert(approxeq(expected, sol));
        }

        template<class MultigridT>
        void test_block_mg(MultigridT &multigrid, const bool verbose = false)
        {
            if(mpi_world_size() > 1) return;

            DVectord rhs;
            DSMatrixd A, I;

            const std::string data_path = Utopia::instance().get("data_path");
            const std::string folder = data_path + "/mg";

            read(folder + "/rhs.bin", rhs);
            read(folder + "/A.bin", A);
            read(folder + "/I.bin", I);

            std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;
            interpolation_operators.push_back(make_ref(I));
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(20);
            multigrid.atol(1e-15);
            multigrid.stol(1e-15);
            multigrid.rtol(1e-15);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(verbose);
            multigrid.fix_semidefinite_operators(true);

            DVectord x = zeros(A.size().get(0));

            int block_size = 2;
            multigrid.block_size(block_size);
            multigrid.update(make_ref(A));

            if(verbose) {
                multigrid.describe();
            }

            multigrid.apply(rhs, x);

            utopia_test_assert(approxeq(rhs, A * x, 1e-6));
        }

        void petsc_block_mg_exp()
        {
           Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> multigrid;
           test_block_mg(multigrid, false);
        }

        void petsc_block_mg()
        {
            Multigrid<DSMatrixd, DVectord> multigrid(
                std::make_shared<SOR<DSMatrixd, DVectord>>(),
                std::make_shared<Factorization<DSMatrixd, DVectord>>()
            );

           // multigrid.set_use_line_search(true);
            // multigrid.set_fix_semidefinite_operators(true);
           test_block_mg(multigrid, false);
        }

        void petsc_cg_mg()
        {
            //! [MG solve example]
            const bool verbose = false;

            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            //reading data from disk
            const std::string data_path = Utopia::instance().get("data_path");
            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;

            //interpolation operators from coarse to fine
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            //choose solver for coarse level solution
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();

#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(Solver::mumps(), Solver::lu_decomposition());
#endif //PETSC_HAVE_MUMPS

            //choose smoother
            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            // auto smoother = std::make_shared<PointJacobi<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(std::move(interpolation_operators));
            multigrid.fix_semidefinite_operators(true);
            multigrid.must_generate_masks(true);
            multigrid.max_it(1);
            multigrid.mg_type(1);
            multigrid.verbose(verbose);
            // multigrid.set_use_line_search(true);

            ConjugateGradient<DSMatrixd, DVectord, HOMEMADE> cg; //with the HOMEMADE works in parallel
            // ConjugateGradient<DSMatrixd, DVectord> cg; //FIXME does not work with the KSP implementation
            cg.verbose(verbose);

            DVectord x_0 = zeros(A.size().get(0));

            //CG with diagonal preconditioner
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<DSMatrixd, DVectord> >());

            x_0 = rhs;
            cg.solve(A, rhs, x_0);

            //CG with multigrid preconditioner
            x_0 = zeros(A.size().get(0));
            cg.set_preconditioner(make_ref(multigrid));
            cg.atol(1e-18);
            cg.rtol(1e-18);
            cg.stol(1e-18);

            x_0 = rhs;
            cg.solve(A, rhs, x_0);

            utopia_test_assert( approxeq(A*x_0, rhs, 1e-6) );

            //Multigrid only
            // x_0 = zeros(A.size().get(0));
            // multigrid.max_it(12);
            // multigrid.verbose(verbose);
            // multigrid.solve(rhs, x_0);

            //! [MG solve example]
        }

        void petsc_mg_jacobi()
        {
            const std::string data_path = Utopia::instance().get("data_path");
            DSMatrixd A, I_1, I_2, I_3;
            DVectord rhs;

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            auto s_A     = local_size(A);
            DVectord x   = local_zeros(s_A.get(0));

            std::vector<std::shared_ptr <DSMatrixd> > interpolation_operators;
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
            auto smoother      = std::make_shared<PointJacobi<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.update(make_ref(A));

            // multigrid.verbose(true);
            // multigrid.set_use_line_search(true);
            multigrid.apply(rhs, x);

            utopia_test_assert( approxeq(A*x, rhs, 1e-6) );
        }

        void petsc_superlu_mg()
        {
            const bool verbose = false;
            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            if(verbose) {
                disp("n unknowns:");
                disp(size(rhs));
            }

            std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;

            // from coarse to fine
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            //  init

#ifdef PETSC_HAVE_SUPERLU_DIST
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >(MATSOLVERSUPERLU_DIST, PCLU);
#else
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >(MATSOLVERPETSC, PCLU);

            if(mpi_world_size() > 1) {
                if(mpi_world_rank() == 0) {
                    std::cerr << "[Error] Direct solver does not work in parallel compile with SuperLU" << std::endl;
                }
                return;
            }
#endif //PETSC_HAVE_SUPERLU_DIST

            // direct_solver->describe(std::cout);

            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(std::move(interpolation_operators));
            multigrid.update(make_ref(A));

            multigrid.max_it(1);
            multigrid.mg_type(2);

            DVectord x_0;

            //! [KSPSolver solve example1]
            KSPSolver<DSMatrixd, DVectord> utopia_ksp;
            auto precond = std::make_shared< InvDiagPreconditioner<DSMatrixd, DVectord> >();
            // utopia_ksp.set_preconditioner(precond);
            utopia_ksp.verbose(verbose);
            utopia_ksp.ksp_type("gmres");
            utopia_ksp.solve(A, rhs, x_0);
            utopia_test_assert( approxeq(A*x_0, rhs, 1e-6) );
            //! [KSPSolver solve example1]

            x_0 = zeros(A.size().get(0));
            multigrid.verbose(false);
            multigrid.max_it(1);
            multigrid.mg_type(1);
            multigrid.pre_smoothing_steps(1);
            multigrid.post_smoothing_steps(1);
            utopia_ksp.set_preconditioner(make_ref(multigrid));
            // utopia_ksp.verbose(true);

            utopia_ksp.atol(1e-18);
            utopia_ksp.rtol(1e-18);
            utopia_ksp.stol(1e-16);

            utopia_ksp.solve(A, rhs, x_0);

            double diff = norm2(rhs - A * x_0);

            if(diff > 1e-6) {
                utopia_error("petsc_superlu_mg fails. Known problem that needs to be fixed!");
            }
            // utopia_test_assert( diff < 1e-6 );
        }


        void petsc_factorization()
        {
            if(mpi_world_size() > 1)
                return;

            DVectord rhs, x;
            DSMatrixd A = zeros(_n, _n);

            assemble_laplacian_1D(A, true);

            x 	= local_zeros(local_size(A).get(0));
            rhs = local_values(local_size(A).get(0), 13.0);

            auto cholesky_factorization = std::make_shared<Factorization<DSMatrixd, DVectord> >();
            cholesky_factorization->set_type(Solver::petsc(), Solver::lu_decomposition());

            if(!cholesky_factorization->solve(A, rhs, x)) {
                utopia_test_assert(false && "failed to solve");
            }

            double diff = norm2(rhs - A * x);
            // utopia_test_assert( approxeq(A * x, rhs, 1e-6) );

            if(diff > 1e-6) {
                utopia_error("petsc_factorization fails. Known problem that needs to be fixed!");
            }
        }


        void petsc_st_cg_mg()
        {
            //! [MG solve example]
            const bool verbose = false;

            DVectord rhs;
            DSMatrixd A, I_1, I_2, I_3;

            //reading data from disk
            const std::string data_path = Utopia::instance().get("data_path");
            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;

            //interpolation operators from coarse to fine
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            //choose solver for coarse level solution
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();

#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(Solver::mumps(), Solver::lu_decomposition());
#endif //PETSC_HAVE_MUMPS

            //choose smoother
            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            // auto smoother = std::make_shared<PointJacobi<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(std::move(interpolation_operators));
            multigrid.max_it(1);
            multigrid.mg_type(1);
            multigrid.verbose(verbose);
            // multigrid.set_use_line_search(true);

            SteihaugToint<DSMatrixd, DVectord, HOMEMADE> cg;
            cg.verbose(verbose);

            DVectord x_0 = zeros(A.size().get(0));

            // plain cg
            cg.atol(1e-18);
            cg.rtol(1e-18);
            cg.stol(1e-18);
            cg.solve(A, -1.0 * rhs, x_0);


            //CG with diagonal preconditioner
            x_0 = zeros(A.size().get(0));
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<DSMatrixd, DVectord> >());
            cg.solve(A, -1.0 * rhs, x_0);


            //CG with multigrid preconditioner
            x_0 = zeros(A.size().get(0));
            cg.set_preconditioner(make_ref(multigrid));
            cg.solve(A, rhs, x_0);


            utopia_test_assert( approxeq(A*x_0, rhs, 1e-6) );

            //Multigrid only
            // x_0 = zeros(A.size().get(0));
            // multigrid.max_it(12);
            // multigrid.verbose(verbose);
            // multigrid.solve(rhs, x_0);

            //! [MG solve example]
        }


        PetscLinearSolverTest()
        : _n(10) { }

    private:
        int _n;
    };

#endif //WITH_PETSC

    void runPetscLinearSolversTest()
    {
        UTOPIA_UNIT_TEST_BEGIN("PetscLinearSolverTest");
#ifdef WITH_PETSC
        PetscLinearSolverTest().run();
#endif

        UTOPIA_UNIT_TEST_END("PetscLinearSolverTest");
    }
}
