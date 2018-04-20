/*
 * @Author: kopanicakova
 * @Date:   2018-02-06 17:47:26
 * @Last Modified by:   kopanicakova
 * @Last Modified time: 2018-04-09 14:03:17
 */
#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

namespace utopia
{

#ifdef WITH_PETSC
    class PetscLinearSolverTest {
    public:
        
        void run()
        {
            UTOPIA_RUN_TEST(petsc_mg_exp);
            UTOPIA_RUN_TEST(petsc_bicgstab);
            UTOPIA_RUN_TEST(petsc_gmres);
            UTOPIA_RUN_TEST(petsc_mg);
            UTOPIA_RUN_TEST(petsc_cg_mg);
            UTOPIA_RUN_TEST(petsc_superlu_cg_mg);
            UTOPIA_RUN_TEST(petsc_mg_jacobi);
            UTOPIA_RUN_TEST(petsc_factorization);
            UTOPIA_RUN_TEST(petsc_block_mg_exp);
            UTOPIA_RUN_TEST(petsc_block_mg);
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
            assert(err < 1e-6);
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
            direct_solver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
#endif //PETSC_HAVE_MUMPS
            
            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            // multigrid.set_use_line_search(true);
            
            
            multigrid.set_transfer_operators(std::move(interpolation_operators));
            // multigrid.set_fix_semidefinite_operators(true);
            multigrid.update(make_ref(A));
            
            DVectord x_0 = zeros(A.size().get(0));
            
            Parameters params;
            params.linear_solver_verbose(false);
            multigrid.set_parameters(params);
            
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
            // gmres->verbose(true);
            gmres->solve(A, rhs, x_0);
            
            assert( approxeq(A*x_0, rhs, 1e-6) );
        }
        
        void petsc_bicgstab()
        {
            
            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);
            
            BiCGStab<DMatrixd, DVectord> bicgs;
            bicgs.solve(mat, rhs, sol);
            
            DVectord expected = zeros(_n);
            assert(approxeq(expected, sol));
        }
        
        void petsc_gmres()
        {
            
            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);
            
            GMRES<DMatrixd, DVectord> gmres;
            gmres.solve(mat, rhs, sol);
            
            DVectord expected = zeros(_n);
            assert(approxeq(expected, sol));
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
            
            multigrid.set_transfer_operators(std::move(interpolation_operators));
            multigrid.max_it(20);
            multigrid.atol(1e-15);
            multigrid.stol(1e-15);
            multigrid.rtol(1e-15);
            multigrid.verbose(verbose);
            
            DVectord x = zeros(A.size().get(0));

            int block_size = 1;
            multigrid.block_size(block_size);
            multigrid.update(make_ref(A));
           
            if(verbose) {
                multigrid.describe();
            }

            multigrid.apply(rhs, x);

            assert(approxeq(rhs, A * x, 1e-6));
        }

        void petsc_block_mg_exp()
        {
           Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> multigrid;
           test_block_mg(multigrid, false);
        }

        void petsc_block_mg()
        {
            Multigrid<DSMatrixd, DVectord> multigrid(
                // std::make_shared<GMRES<DSMatrixd, DVectord>>(),
                std::make_shared<GaussSeidel<DSMatrixd, DVectord>>(),
                // std::make_shared<Factorization<DSMatrixd, DVectord>>()
                std::make_shared<LUDecomposition<DSMatrixd, DVectord>>()
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
            direct_solver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
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
            
            // ConjugateGradient<DSMatrixd, DVectord, HOMEMADE> cg; //with the HOMEMADE works in parallel
            ConjugateGradient<DSMatrixd, DVectord> cg;
            cg.verbose(verbose);
            
            DVectord x_0 = zeros(A.size().get(0));
            
            //CG with diagonal preconditioner
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<DSMatrixd, DVectord> >());
            cg.solve(A, rhs, x_0);
            
            //CG with multigrid preconditioner
            x_0 = zeros(A.size().get(0));
            cg.set_preconditioner(make_ref(multigrid));
            cg.atol(1e-18);
            cg.rtol(1e-18);
            cg.stol(1e-18);
            
            cg.solve(A, rhs, x_0);
            
            assert( approxeq(A*x_0, rhs, 1e-6) );
            
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
            
            assert( approxeq(A*x, rhs, 1e-6) );
        }
        
        void petsc_superlu_cg_mg()
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
            auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
#ifdef PETSC_HAVE_SUPERLU_DIST
            direct_solver->set_type(SUPERLU_DIST_TAG, LU_DECOMPOSITION_TAG);
#else
            if(mpi_world_size() > 1) {
                if(mpi_world_rank() == 0) {
                    std::cerr << "[Error] Direct solver does not work in parallel compile with SuperLU" << std::endl;
                }
                return;
            }
#endif //PETSC_HAVE_SUPERLU_DIST
            
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
            assert( approxeq(A*x_0, rhs, 1e-6) );
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
                utopia_error("petsc_superlu_cg_mg fails. Known problem that needs to be fixed!");
            }
            // assert( diff < 1e-6 );
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
            cholesky_factorization->set_type(PETSC_TAG, LU_DECOMPOSITION_TAG);
            
            if(!cholesky_factorization->solve(A, rhs, x)) {
                assert(false && "failed to solve");
            }
            
            double diff = norm2(rhs - A * x);
            // assert( approxeq(A * x, rhs, 1e-6) );
            
            if(diff > 1e-6) {
                utopia_error("petsc_factorization fails. Known problem that needs to be fixed!");
            }
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
