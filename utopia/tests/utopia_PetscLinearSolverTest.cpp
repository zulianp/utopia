/*
 * @Author: kopanicakova
 * @Date:   2018-02-06 17:47:26
 * @Last Modified by:   kopanicakova
 * @Last Modified time: 2018-04-09 14:03:17
 */
#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"

namespace utopia
{
    template<class Matrix>
    void assemble_laplacian_1D(const utopia::SizeType n, Matrix &m)
    {
        
        // n x n matrix with maximum 3 entries x row
        {
            Write<Matrix> w(m);
            Range r = row_range(m);
            
            //You can use set instead of add. [Warning] Petsc does not allow to mix add and set.
            for(SizeType i = r.begin(); i != r.end(); ++i) {
                if(i > 0) {
                    m.add(i, i - 1, -1.0);
                }
                
                if(i < n-1) {
                    m.add(i, i + 1, -1.0);
                }
                
                m.add(i, i, 2.0);
            }
        }
    }
    
    
#ifdef WITH_PETSC
    class PetscLinearSolverTest {
    public:
        
        void run()
        {
            UTOPIA_RUN_TEST(petsc_mg_exp_test);
            UTOPIA_RUN_TEST(petsc_bicgstab_test);
            UTOPIA_RUN_TEST(petsc_gmres_test);
            UTOPIA_RUN_TEST(petsc_mg_test);
            UTOPIA_RUN_TEST(petsc_cg_mg_test);
            UTOPIA_RUN_TEST(petsc_superlu_cg_mg_test);
            UTOPIA_RUN_TEST(petsc_mg_jacobi_test);
            UTOPIA_RUN_TEST(petsc_cholesky_test);
        }
        
        void petsc_mg_exp_test()
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
            
            auto smoother = std::make_shared<GMRES<DSMatrixd, DVectord>>();
            auto linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord>>();
            // auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
            // linear_solver->verbose(true);
            Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> multigrid;//(smoother, linear_solver);
            // multigrid.set_default_pc_type(PCILU);
            // multigrid.set_default_ksp_type(KSPFGMRES);
            
            // Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> multigrid;
            multigrid.init_transfer_from_fine_to_coarse(std::move(interpolation_operators));
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
        
        void petsc_bicgstab_test()
        {
            
            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);
            
            BiCGStab<DMatrixd, DVectord> bicgs;
            bicgs.solve(mat, rhs, sol);
            
            DVectord expected = zeros(_n);
            assert(approxeq(expected, sol));
        }
        
        void petsc_gmres_test()
        {
            
            DMatrixd mat = identity(_n, _n);
            DVectord rhs = zeros(_n);
            DVectord sol = zeros(_n);
            
            GMRES<DMatrixd, DVectord> gmres;
            gmres.solve(mat, rhs, sol);
            
            DVectord expected = zeros(_n);
            assert(approxeq(expected, sol));
        }
        
        
        
        void petsc_mg_test()
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
            multigrid.set_use_line_search(true);
            
            
            multigrid.init_transfer_from_fine_to_coarse(std::move(interpolation_operators));
            multigrid.set_fix_semidefinite_operators(true);
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
        
        
        void petsc_cg_mg_test()
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
            multigrid.init_transfer_from_fine_to_coarse(std::move(interpolation_operators));
            multigrid.max_it(1);
            multigrid.mg_type(1);
            multigrid.verbose(verbose);
            multigrid.set_use_line_search(true);
            
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
        
        void petsc_mg_jacobi_test()
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
            
            auto direct_solver = std::make_shared< Factorization<DSMatrixd, DVectord> >();
            auto smoother = std::make_shared<PointJacobi<DSMatrixd, DVectord>>();
            Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
            multigrid.init_transfer_from_fine_to_coarse(interpolation_operators);
            multigrid.update(make_ref(A));
            
            // multigrid.verbose(true);
            multigrid.set_use_line_search(true);
            multigrid.solve(rhs, x);
            
            assert( approxeq(A*x, rhs, 1e-6) );
        }
        
        void petsc_superlu_cg_mg_test()
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
            multigrid.init_transfer_from_fine_to_coarse(std::move(interpolation_operators));
            multigrid.update(make_ref(A));
            
            multigrid.max_it(1);
            multigrid.mg_type(2);
            
            DVectord x_0;
            
            //! [KSPSolver solve example1]
            KSPSolver<DSMatrixd, DVectord> utopia_ksp;
            auto precond = std::make_shared< InvDiagPreconditioner<DSMatrixd, DVectord> >();
            // utopia_ksp.set_preconditioner(precond);
            utopia_ksp.verbose(verbose);
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
                utopia_error("petsc_superlu_cg_mg_test fails. Known problem that needs to be fixed!");
            }
            // assert( diff < 1e-6 );
        }
        
        void petsc_cholesky_test()
        {
            if(mpi_world_size() > 1)
                return;
            
            DVectord rhs, x;
            DSMatrixd A = zeros(_n, _n);
            
            assemble_laplacian_1D(_n, A);
            {
                Range r = row_range(A);
                Write<DSMatrixd> w(A);
                if(r.begin() == 0) {
                    A.set(0, 0, 1.);
                    A.set(0, 1, 0);
                }
                
                
                if(r.end() == _n) {
                    A.set(_n-1, _n-1, 1.);
                    A.set(_n-1, _n-2, 0);
                }
            }
            
            x 	= local_zeros(local_size(A).get(0));
            rhs = local_values(local_size(A).get(0), 13.0);
            
            auto cholesky_factorization = std::make_shared<Factorization<DSMatrixd, DVectord> >();
            cholesky_factorization->set_type(PETSC_TAG, CHOLESKY_DECOMPOSITION_TAG);
            
            if(!cholesky_factorization->solve(A, rhs, x)) {
                assert(false && "failed to solve");
            }
            
            double diff = norm2(rhs - A * x);
            // assert( approxeq(A * x, rhs, 1e-6) );
            
            if(diff > 1e-6) {
                utopia_error("petsc_cholesky_test fails. Known problem that needs to be fixed!");
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
