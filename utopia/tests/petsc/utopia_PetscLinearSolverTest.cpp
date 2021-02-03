#include "utopia.hpp"
#include "utopia_TestProblems.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

namespace utopia {

#ifdef UTOPIA_WITH_PETSC
    class PetscLinearSolverTest {
    public:
        using Traits = utopia::Traits<PetscVector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Comm comm_;

        void run() {
            UTOPIA_RUN_TEST(petsc_cg);
            UTOPIA_RUN_TEST(petsc_bicgstab);
            UTOPIA_RUN_TEST(petsc_gmres);
// FIXME make it work also without mumps
#ifdef PETSC_HAVE_MUMPS
            UTOPIA_RUN_TEST(petsc_mg);
            UTOPIA_RUN_TEST(petsc_cg_mg);
            UTOPIA_RUN_TEST(petsc_mg_1D);

            UTOPIA_RUN_TEST(petsc_block_mg_exp);  // petsc 3.11.3 ERROR here
            UTOPIA_RUN_TEST(petsc_block_mg);
            UTOPIA_RUN_TEST(petsc_mg_exp);
            UTOPIA_RUN_TEST(petsc_superlu_mg);
            UTOPIA_RUN_TEST(petsc_mg_jacobi);
            UTOPIA_RUN_TEST(petsc_factorization);
            UTOPIA_RUN_TEST(petsc_st_cg_mg);
            UTOPIA_RUN_TEST(petsc_redundant_test);

#endif  // PETSC_HAVE_MUMPS
        }

        void petsc_cg() {
            Poisson1D<PetscMatrix, PetscVector> fun(100, 2);
            PetscVector x = fun.initial_guess();
            PetscVector rhs;
            fun.get_rhs(rhs);
            PetscMatrix A;
            fun.hessian(x, A);
            rhs *= 0.00001;

            ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> cg;
            cg.rtol(1e-6);
            cg.atol(1e-6);
            cg.max_it(500);
            // cg.verbose(true);

            x = rhs;
            cg.update(std::make_shared<PetscMatrix>(A));
            cg.apply(rhs, x);

            utopia_test_assert(approxeq(rhs, A * x, 1e-5));
        }

        void petsc_mg_1D() {
            if (comm_.size() > 1) {
                return;
            }

            const static bool verbose = false;

            // MultiLevelTestProblem<PetscMatrix, PetscVector> ml_problem(4, 2);
            // ml_problem.write_matlab("./");

            MultiLevelTestProblem1D<PetscMatrix, PetscVector, Poisson1D<PetscMatrix, PetscVector>> ml_problem(
                4, 10, true);

            auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
            smoother->verbose(verbose);

            Multigrid<PetscMatrix, PetscVector> multigrid(
                smoother, std::make_shared<Factorization<PetscMatrix, PetscVector>>()
                // std::make_shared<ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE>>(),
                // std::make_shared<ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE>>()
            );

            multigrid.set_transfer_operators(ml_problem.get_transfer());
            multigrid.max_it(4);
            multigrid.atol(1e-12);
            multigrid.stol(1e-10);
            multigrid.rtol(1e-10);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(false);

            auto fun = ml_problem.get_functions().back();
            PetscMatrix A;
            PetscVector g, x;
            fun->get_eq_constrains_values(x);
            fun->gradient(x, g);
            fun->hessian(x, A);

            // PetscVector x = zeros(size(*ml_problem.rhs));
            multigrid.update(std::make_shared<PetscMatrix>(A));

            if (verbose) {
                multigrid.describe();
            }

            multigrid.apply(g, x);
            utopia_test_assert(approxeq(g, A * x, 1e-7));
        }

        void petsc_mg_exp() {
            PetscVector rhs;
            PetscMatrix A, I_1, I_2, I_3;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
            auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
            // auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
            // linear_solver->verbose(true);
            Multigrid<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL> multigrid;  //(smoother, linear_solver);
            // multigrid.set_default_pc_type(PCILU);
            // multigrid.set_default_ksp_type(KSPFGMRES);

            // Multigrid<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL> multigrid;
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(20);
            multigrid.atol(1e-15);
            multigrid.stol(1e-15);
            multigrid.rtol(1e-15);

            PetscVector x(row_layout(A), 0.0);
            // multigrid.verbose(true);
            multigrid.solve(A, rhs, x);

            const double err = norm2(A * x - rhs);
            utopia_test_assert(err < 1e-6);
        }

        void petsc_mg() {
            // reading data from outside
            PetscVector rhs;
            PetscMatrix A, I_1, I_2, I_3;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            // read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;

            // from coarse to fine
            // interpolation_operators.push_back(std::move(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            //  init
            auto direct_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(Solver::mumps(), Solver::lu_decomposition());
#endif  // PETSC_HAVE_MUMPS

            auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();

            // auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector, HOMEMADE>>(); smoother->l1(true);

            Multigrid<PetscMatrix, PetscVector> multigrid(smoother, direct_solver);
            // multigrid.set_use_line_search(true);
            // multigrid.verbose(true);

            multigrid.set_transfer_operators(interpolation_operators);
            // multigrid.set_fix_semidefinite_operators(true);
            multigrid.update(make_ref(A));

            PetscVector x_0(row_layout(A), 0.0);

            // multigrid.verbose(true);
            multigrid.apply(rhs, x_0);

            double diff = norm2(A * x_0 - rhs);
            utopia_test_assert(diff < 1e-6);
            // std::cout<<"diff: "<< diff << " \n";

            // multigrid.verbose(false);

            x_0.zeros(row_layout(A));
            multigrid.cycle_type(FULL_CYCLE);
            multigrid.apply(rhs, x_0);

            x_0.zeros(row_layout(A));
            multigrid.cycle_type(FULL_CYCLE);
            multigrid.v_cycle_repetition(2);

            multigrid.apply(rhs, x_0);
            diff = norm2(A * x_0 - rhs);
            utopia::out() << "diff: " << diff << " \n";

            multigrid.max_it(1);
            multigrid.cycle_type(MULTIPLICATIVE_CYCLE);
            auto gmres = std::make_shared<GMRES<PetscMatrix, PetscVector>>();
            gmres->set_preconditioner(make_ref(multigrid));

            x_0 = 0.0 * x_0;
            gmres->verbose(false);
            gmres->atol(1e-16);
            gmres->rtol(1e-16);
            gmres->solve(A, rhs, x_0);

            // std::cout<<"gmres.getn_um_it(): "<< gmres->getn_um_it() << "   \n";

            if (diff > 1e-6) {
                utopia_error("petsc_mg: gmres preconditioned with mg does not do what it is supposed to");
            }

            // utopia_test_assert( approxeq(diff, 0., 1e-6) );
        }

        void petsc_bicgstab() {
            PetscMatrix mat;
            mat.identity(layout(comm_, Traits::decide(), Traits::decide(), n_, n_), 1.0);

            PetscVector rhs(row_layout(mat), 0.0);
            PetscVector sol(row_layout(mat), 0.0);

            BiCGStab<PetscMatrix, PetscVector> bicgs;
            bicgs.solve(mat, rhs, sol);

            PetscVector expected(row_layout(mat), 0.0);
            utopia_test_assert(approxeq(expected, sol));
        }

        void petsc_gmres() {
            PetscMatrix mat;
            mat.identity(layout(comm_, Traits::decide(), Traits::decide(), n_, n_), 1.0);

            PetscVector rhs(row_layout(mat), 0.0);
            PetscVector sol(row_layout(mat), 0.0);

            GMRES<PetscMatrix, PetscVector> gmres;
            gmres.pc_type(PCBJACOBI);

            gmres.number_of_subdomains(comm_.size());
            gmres.update(std::make_shared<PetscMatrix>(mat));
            gmres.sub_ksp_pc_type(KSPPREONLY, PCILU);

            gmres.verbose(false);
            gmres.apply(rhs, sol);

            PetscVector expected(row_layout(mat), 0.0);
            utopia_test_assert(approxeq(expected, sol));
        }

        void petsc_redundant_test() {
            PetscMatrix mat;
            mat.identity(layout(comm_, Traits::decide(), Traits::decide(), n_, n_), 1.0);

            PetscVector rhs(row_layout(mat), 0.0);
            PetscVector sol(row_layout(mat), 0.0);

            auto solver = std::make_shared<utopia::RedundantLinearSolver<PetscMatrix, PetscVector>>();
            solver->number_of_parallel_solves(comm_.size());
            // solver->number_of_parallel_solves(1);
            solver->ksp_type("gmres");
            solver->pc_type("lu");

            solver->verbose(false);
            solver->solve(mat, rhs, sol);

            PetscVector expected(row_layout(mat), 0.0);
            utopia_test_assert(approxeq(expected, sol));
        }

        template <class MultigridT>
        void test_block_mg(MultigridT &multigrid, const bool verbose) {
            // FiXME USE serial comm
            if (comm_.size() > 1) {
                return;
            }

            PetscVector rhs;
            PetscMatrix A, I;

            const std::string data_path = Utopia::instance().get("data_path");
            const std::string folder = data_path + "/mg";

            read(folder + "/rhs.bin", rhs);
            read(folder + "/A.bin", A);
            read(folder + "/I.bin", I);

            std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;
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

            PetscVector x(row_layout(A), 0.0);

            int block_size = 2;
            multigrid.block_size(block_size);
            multigrid.update(make_ref(A));

            if (verbose) {
                multigrid.describe();
            }

            multigrid.apply(rhs, x);

            utopia_test_assert(approxeq(rhs, A * x, 1e-6));
        }

        void petsc_block_mg_exp() {
            Multigrid<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL> multigrid;
            test_block_mg(multigrid, false);
        }

        void petsc_block_mg() {
            Multigrid<PetscMatrix, PetscVector> multigrid(std::make_shared<SOR<PetscMatrix, PetscVector>>(),
                                                          std::make_shared<Factorization<PetscMatrix, PetscVector>>());

            // multigrid.set_use_line_search(true);
            // multigrid.set_fix_semidefinite_operators(true);
            test_block_mg(multigrid, false);
        }

        void petsc_cg_mg() {
            //! [MG solve example]
            const bool verbose = false;

            PetscVector rhs;
            PetscMatrix A, I_1, I_2, I_3;

            // reading data from disk
            const std::string data_path = Utopia::instance().get("data_path");
            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;

            // interpolation operators from coarse to fine
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            // choose solver for coarse level solution
            auto direct_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();

#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(Solver::mumps(), Solver::lu_decomposition());
#endif  // PETSC_HAVE_MUMPS

            // choose smoother
            auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
            // auto smoother = std::make_shared<PointJacobi<PetscMatrix, PetscVector>>();
            Multigrid<PetscMatrix, PetscVector> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.fix_semidefinite_operators(true);
            multigrid.must_generate_masks(true);
            multigrid.max_it(1);
            multigrid.mg_type(true);
            multigrid.verbose(verbose);
            // multigrid.set_use_line_search(true);

            ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> cg;  // with the HOMEMADE works in parallel
            // ConjugateGradient<PetscMatrix, PetscVector> cg; //FIXME does not work with the KSP implementation
            cg.verbose(verbose);

            PetscVector x_0(row_layout(A), 0.0);

            // CG with diagonal preconditioner
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>());

            x_0 = rhs;
            cg.solve(A, rhs, x_0);

            // CG with multigrid preconditioner
            x_0.zeros(row_layout(A));
            cg.set_preconditioner(make_ref(multigrid));
            cg.atol(1e-18);
            cg.rtol(1e-18);
            cg.stol(1e-18);

            x_0 = rhs;
            cg.solve(A, rhs, x_0);

            if (!approxeq(A * x_0, rhs, 1e-6)) {
                m_utopia_error("utopia_test_assert( approxeq(A*x_0, rhs, 1e-6) )");
            }
            // utopia_test_assert( approxeq(A*x_0, rhs, 1e-6) );

            // Multigrid only
            // PetscVector x_0(row_layout(A), 0.0);
            // multigrid.max_it(12);
            // multigrid.verbose(verbose);
            // multigrid.solve(rhs, x_0);

            //! [MG solve example]
        }

        void petsc_mg_jacobi() {
            const std::string data_path = Utopia::instance().get("data_path");
            PetscMatrix A, I_1, I_2, I_3;
            PetscVector rhs;

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            PetscVector x(row_layout(A), 0.0);

            std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            auto direct_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
            auto smoother = std::make_shared<PointJacobi<PetscMatrix, PetscVector>>();
            Multigrid<PetscMatrix, PetscVector> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.update(make_ref(A));

            // multigrid.verbose(true);
            // multigrid.set_use_line_search(true);
            multigrid.apply(rhs, x);

            utopia_test_assert(approxeq(A * x, rhs, 1e-6));
        }

        void petsc_superlu_mg() {
            const bool verbose = false;
            PetscVector rhs;
            PetscMatrix A, I_1, I_2, I_3;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            if (verbose) {
                disp("n unknowns:");
                disp(size(rhs));
            }

            std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;

            // from coarse to fine
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            //  init

#ifdef PETSC_HAVE_SUPERLU_DIST
            auto direct_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>(MATSOLVERSUPERLU_DIST, PCLU);
#else
            auto direct_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>(MATSOLVERPETSC, PCLU);

            if (comm_.size() > 1) {
                if (mpi_world_rank() == 0) {
                    std::cerr << "[Error] Direct solver does not work in parallel compile with SuperLU" << std::endl;
                }
                return;
            }
#endif  // PETSC_HAVE_SUPERLU_DIST

            // direct_solver->describe(std::cout);

            auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
            Multigrid<PetscMatrix, PetscVector> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.update(make_ref(A));

            multigrid.max_it(1);
            multigrid.mg_type(true);

            PetscVector x_0(row_layout(A), 0.0);

            //! [KSPSolver solve example1]
            KSPSolver<PetscMatrix, PetscVector> utopia_ksp;
            auto precond = std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>();
            // utopia_ksp.set_preconditioner(precond);
            utopia_ksp.verbose(verbose);
            utopia_ksp.ksp_type("gmres");
            utopia_ksp.solve(A, rhs, x_0);
            utopia_test_assert(approxeq(A * x_0, rhs, 1e-6));
            //! [KSPSolver solve example1]

            multigrid.verbose(verbose);
            multigrid.max_it(1);
            multigrid.mg_type(true);
            multigrid.pre_smoothing_steps(1);
            multigrid.post_smoothing_steps(1);
            utopia_ksp.set_preconditioner(make_ref(multigrid));
            utopia_ksp.verbose(verbose);

            utopia_ksp.atol(1e-18);
            utopia_ksp.rtol(1e-18);
            utopia_ksp.stol(1e-16);
            utopia_ksp.max_it(10);
            utopia_ksp.norm_type("preconditioned");

            utopia_ksp.solve(A, rhs, x_0);

            double diff = norm2(A * x_0 - rhs);

            if (diff > 1e-6) {
                utopia_error("petsc_superlu_mg fails. Known problem that needs to be fixed!");
            }

            utopia_test_assert(diff < 1e-6);
        }

        void petsc_factorization() {
            if (comm_.size() > 1) {
                return;
            }

            PetscVector rhs, x;
            PetscMatrix A;
            A.dense(serial_layout(n_, n_), 0.0);

            assemble_laplacian_1D(A, true);

            x.zeros(row_layout(A));
            rhs.values(row_layout(A), 13.0);

            auto cholesky_factorization = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
            cholesky_factorization->set_type(Solver::petsc(), Solver::lu_decomposition());

            if (cholesky_factorization->solve(A, rhs, x)) {
                double diff = norm2(rhs - A * x);
                if (diff > 1e-6) {
                    utopia_error("petsc_factorization fails. Known problem that needs to be fixed!");
                }

            } else {
                utopia_error("petsc_factorization fails. Known problem that needs to be fixed!");
            }
        }

        void petsc_st_cg_mg() {
            //! [MG solve example]
            const bool verbose = false;

            PetscVector rhs;
            PetscMatrix A, I_1, I_2, I_3;

            // reading data from disk
            const std::string data_path = Utopia::instance().get("data_path");
            read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
            read(data_path + "/laplace/matrices_for_petsc/f_A", A);
            read(data_path + "/laplace/matrices_for_petsc/I_1", I_1);
            read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
            read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);

            std::vector<std::shared_ptr<PetscMatrix>> interpolation_operators;

            // interpolation operators from coarse to fine
            interpolation_operators.push_back(make_ref(I_1));
            interpolation_operators.push_back(make_ref(I_2));
            interpolation_operators.push_back(make_ref(I_3));

            // choose solver for coarse level solution
            auto direct_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();

#ifdef PETSC_HAVE_MUMPS
            direct_solver->set_type(Solver::mumps(), Solver::lu_decomposition());
#endif  // PETSC_HAVE_MUMPS

            // choose smoother
            auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
            // auto smoother = std::make_shared<PointJacobi<PetscMatrix, PetscVector>>();
            Multigrid<PetscMatrix, PetscVector> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(1);
            multigrid.mg_type(true);
            multigrid.verbose(verbose);
            // multigrid.set_use_line_search(true);

            SteihaugToint<PetscMatrix, PetscVector, HOMEMADE> cg;
            cg.verbose(verbose);
            cg.atol(1e-14);

            PetscVector x_0(row_layout(A), 0.0);

            // plain cg
            cg.atol(1e-13);
            cg.rtol(1e-13);
            cg.stol(1e-13);
            cg.verbose(verbose);
            cg.solve(A, rhs, x_0);

            // CG with diagonal preconditioner
            x_0.set(0.0);
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>());
            cg.solve(A, rhs, x_0);

            // CG with multigrid preconditioner
            x_0.set(0.0);
            cg.set_preconditioner(make_ref(multigrid));
            cg.max_it(50);
            cg.solve(A, rhs, x_0);

            utopia_test_assert(approxeq(A * x_0, rhs, 1e-6));

            // Multigrid only
            // PetscVector x_0(row_layout(A), 0.0);
            // multigrid.max_it(12);
            // multigrid.verbose(verbose);
            // multigrid.solve(rhs, x_0);

            //! [MG solve example]
        }

        PetscLinearSolverTest() = default;

    private:
        int n_{10};
    };

#endif  // UTOPIA_WITH_PETSC

    static void petsc_linear() {
#ifdef UTOPIA_WITH_PETSC
        PetscLinearSolverTest().run();
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(petsc_linear);
}  // namespace utopia
