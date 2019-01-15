#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "utopia_InputParameters.hpp"

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



#ifdef  WITH_SLEPC
	class SlepcsSolverTest {
	public:
		
		void run()
		{
			UTOPIA_RUN_TEST(petsc_slepc_generalized_eigen_test); 
			UTOPIA_RUN_TEST(petsc_slepc_eigen_test); 
			UTOPIA_RUN_TEST(tr_more_sorensen_eigen_test); 
			UTOPIA_RUN_TEST(pseudo_tr_test); 
			UTOPIA_RUN_TEST(pseudo_tr_stiff_test); 
			UTOPIA_RUN_TEST(pseudo_cont_test); 
			UTOPIA_RUN_TEST(lm_test); 
			UTOPIA_RUN_TEST(rosenbrock_test); 
		}



		void petsc_slepc_eigen_test()
		{

			DSMatrixd A = sparse(_n, _n, 3);
			assemble_laplacian_1D(_n, A);

			bool verbose = false; 

			SlepcSolver<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> slepc; 

			slepc.portion_of_spectrum("smallest_real"); 
			slepc.verbose(verbose); 
			slepc.solve(A); 
			slepc.print_eigenpairs(); 

			PetscScalar iegr, eigi; 
			DVectord vr, vi; 

			slepc.get_eigenpairs(0, iegr, eigi, vr, vi); 

			slepc.get_real_eigenpair(1, iegr, vr); 

		}

		void petsc_slepc_generalized_eigen_test()
		{
			DSMatrixd A = sparse(_n, _n, 3);
			assemble_laplacian_1D(_n, A);

			DSMatrixd B = 99 * local_identity(local_size(A)); 
			bool verbose = false; 

			SlepcSolver<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> slepc; 

			slepc.portion_of_spectrum("largest_real"); 
			slepc.verbose(verbose); 
			slepc.problem_type("generalized_hermitian");
			slepc.solver_type("krylovschur"); 
			slepc.solve(A, B); 
			slepc.print_eigenpairs(); 

			PetscScalar iegr, eigi; 
			DVectord vr, vi; 

			slepc.get_eigenpairs(0, iegr, eigi, vr, vi); 
			slepc.get_real_eigenpair(1, iegr, vr); 
		}

	    void tr_more_sorensen_eigen_test()
	    {
	    	if(mpi_world_size() != 1)
	    		return; 

	       	DVectord x_w1  = values(4, 10);
			DVectord expected_woods = values(4, 1);
			Woods14<DMatrixd, DVectord> fun_woods;

			bool verbose = false; 

			{
				Write<DVectord> w1(x_w1);
				x_w1.set(0, -3);
				x_w1.set(1, -1);
				x_w1.set(2, -3);
				x_w1.set(3, -1);
			}
					
			auto subproblem = std::make_shared<utopia::KSP_TR<DMatrixd, DVectord> >();
			subproblem->ksp_type("gltr"); 
			subproblem->atol(1e-10); 
			subproblem->rtol(1e-10); 

			TrustRegion<DMatrixd, DVectord> tr_solver(subproblem);
			tr_solver.verbose(verbose);
			tr_solver.max_it(500); 
			tr_solver.atol(1e-11); 
			tr_solver.rtol(1e-14); 
			tr_solver.stol(1e-12); 
			tr_solver.solve(fun_woods, x_w1);				

			auto eigen_solver = std::make_shared<SlepcSolver<DMatrixd, DVectord, PETSC_EXPERIMENTAL> >();
			
			#ifdef SLEPC_HAVE_ARPACK
				eigen_solver->solver_type("arpack");			
			#endif
			
			auto linear_solver = std::make_shared<LUDecomposition<DMatrixd, DVectord> >();
			linear_solver->set_library_type(PETSC_TAG); 

			auto ms_subproblem = std::make_shared<utopia::MoreSorensenEigen<DMatrixd, DVectord> >(linear_solver, eigen_solver);


			// reset initial guess
			{
				Write<DVectord> w1(x_w1);
				x_w1.set(0, -3);
				x_w1.set(1, -1);
				x_w1.set(2, -3);
				x_w1.set(3, -1);
			}

			TrustRegion<DMatrixd, DVectord> tr_solver2(ms_subproblem);
			tr_solver2.verbose(verbose);
			tr_solver2.max_it(500); 
			tr_solver2.atol(1e-11); 
			tr_solver2.rtol(1e-14); 
			tr_solver2.stol(1e-15); 
			tr_solver2.solve(fun_woods, x_w1);				
	    }

		void pseudo_tr_test()
		{
	    	if(mpi_world_size() != 1)
	    		return; 

			DVectord x  = values(4, 10);
			DVectord expected_woods = values(4, 1);
			Woods14<DMatrixd, DVectord> fun;

			auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord>>();	
			linear_solver->atol(1e-14); 
			linear_solver->max_it(10000);

			auto eigen_solver = std::make_shared<SlepcSolver<DMatrixd, DVectord, PETSC_EXPERIMENTAL> >();
			
			#ifdef SLEPC_HAVE_ARPACK
				eigen_solver->solver_type("arpack");			
			#endif		

			PseudoTrustRegion<DMatrixd, DVectord> solver(linear_solver, eigen_solver); 

			solver.atol(1e-9); 
			solver.stol(1e-14); 
			solver.max_it(500);
			solver.verbose(false); 

			solver.solve(fun, x); 
			utopia_test_assert(approxeq(x, expected_woods, 1e-8));
		}


		void pseudo_tr_stiff_test()
		{
			const SizeType n = 100; 

			MildStiffExample<DMatrixd, DVectord> fun(n); 
			DVectord x; 
			fun.get_initial_guess(x); 

			auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord>>();	
			linear_solver->atol(1e-14); 
			linear_solver->max_it(10000);

			auto eigen_solver = std::make_shared<SlepcSolver<DMatrixd, DVectord, PETSC_EXPERIMENTAL> >();
			
			#ifdef SLEPC_HAVE_ARPACK
				eigen_solver->solver_type("arpack");			
			#endif		

			PseudoTrustRegion<DMatrixd, DVectord> solver(linear_solver, eigen_solver); 

			solver.atol(1e-9); 
			solver.stol(1e-14); 
			solver.max_it(1);
			solver.verbose(false); 

			solver.solve(fun, x); 
		}


		void pseudo_cont_test()
		{
			if(mpi_world_size() != 1)
	    		return; 

			DVectord x  = values(4, 10);
			DVectord expected_woods = values(4, 1);
			Woods14<DMatrixd, DVectord> fun;

			auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord>>();	
			linear_solver->atol(1e-14); 
			linear_solver->max_it(10000);

			PseudoContinuation<DMatrixd, DVectord> solver(linear_solver); 

			solver.atol(1e-9); 
			solver.stol(1e-14); 
			solver.max_it(500);
			solver.verbose(false); 

			solver.solve(fun, x); 
			utopia_test_assert(approxeq(x, expected_woods, 1e-8));


			const SizeType n = 100; 

			MildStiffExample<DMatrixd, DVectord> fun_stiff(n); 
			DVectord x_stiff; 
			fun_stiff.get_initial_guess(x_stiff); 
			solver.solve(fun_stiff, x_stiff); 

		}

		void lm_test()
		{
			const SizeType n = 100; 

			MildStiffExample<DMatrixd, DVectord> fun_stiff(n); 
			DVectord x_stiff; 
			fun_stiff.get_initial_guess(x_stiff); 

			auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord>>();	
			linear_solver->atol(1e-14); 
			linear_solver->max_it(10000);
			LevenbergMarquardt<DMatrixd, DVectord> solver(linear_solver); 

			solver.atol(1e-9); 
			solver.stol(1e-14); 
			solver.max_it(500);
			solver.verbose(false); 

			auto params_ls = std::make_shared<InputParameters>(); 
			params_ls->set("atol", 1e-12); 
			auto params_ls_cast = std::static_pointer_cast<Input>(params_ls); 

			InputParameters in;
			in.set("tau0", 1);
			in.set("linear-solver", params_ls_cast);


			solver.read(in);
			// solver.print_usage(std::cout); 
			solver.solve(fun_stiff, x_stiff); 

		}

		void rosenbrock_test()
		{
			const SizeType n = 100; 

			MildStiffExample<DMatrixd, DVectord> fun_stiff(n); 
			DVectord x_stiff; 
			fun_stiff.get_initial_guess(x_stiff); 

			auto linear_solver = std::make_shared<GMRES<DMatrixd, DVectord>>();	
			linear_solver->atol(1e-14); 
			linear_solver->max_it(10000);
			
			auto eigen_solver = std::make_shared<SlepcSolver<DMatrixd, DVectord, PETSC_EXPERIMENTAL> >();
			
			#ifdef SLEPC_HAVE_ARPACK
				eigen_solver->solver_type("arpack");			
			#endif

			RosenbrockTrustRegion<DMatrixd, DVectord> solver(linear_solver, eigen_solver); 

			solver.atol(1e-12); 
			solver.stol(1e-14); 
			solver.max_it(100);
			solver.verbose(false); 

			solver.solve(fun_stiff, x_stiff); 

			if(mpi_world_size() != 1)
	  		  		return; 

			DVectord x  = values(4, 10);
			DVectord expected_woods = values(4, 1);
			Woods14<DMatrixd, DVectord> fun;

			solver.solve(fun, x); 

		}

		SlepcsSolverTest()
		: _n(10) { }
		
	private:
		int _n;
	};

#endif //WITH_SLEPC


	void runPetscSlepcSolversTest()
	{
		UTOPIA_UNIT_TEST_BEGIN("runSlepcsSolverTest");
#ifdef  WITH_SLEPC
				SlepcsSolverTest().run();
#endif		
		UTOPIA_UNIT_TEST_END("runSlepcsSolverTest");				
	}
}
