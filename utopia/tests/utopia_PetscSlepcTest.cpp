/*
* @Author: kopanicakova
* @Date:   2018-02-06 17:47:26
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-04-11 14:55:23
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



#ifdef  WITH_SLEPC
	class SlepcsSolverTest {
	public:
		
		void run()
		{
			UTOPIA_RUN_TEST(petsc_slepc_generalized_eigen_test); 
			UTOPIA_RUN_TEST(petsc_slepc_eigen_test); 
			UTOPIA_RUN_TEST(nested_mat_test); 
			UTOPIA_RUN_TEST(tr_more_sorensen_eigen_test); 
		}



		void petsc_slepc_eigen_test()
		{

			DSMatrixd A = sparse(_n, _n, 3);
			assemble_laplacian_1D(_n, A);

			bool verbose = false; 

			EigenValueSlover<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> slepc; 

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

			EigenValueSlover<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> slepc; 

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


	    void nested_mat_test()
	    {
	        DSMatrixd merged_mat; 
	        DVectord merged_vec; 

	        DSMatrixd M00 = local_identity(20, 20);
	        DSMatrixd M01 = 2* local_identity(20, 30);
	        DSMatrixd M10 = 3* local_identity(30, 20);
	        DSMatrixd M11 = 4 *local_identity(30, 30);

	        Mat matrices[4]; 
	        matrices[0] = raw_type(M00); 
	        matrices[1] = raw_type(M01); 
	        matrices[2] = raw_type(M10); 
	        matrices[3] = raw_type(M11); 

	        MatDestroy(&raw_type(merged_mat)); 
	       	MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, matrices, &raw_type(merged_mat));

	        Vec x[2]; 

	        DVectord x1 = local_values(20, 3.0); 
	        DVectord x2 = local_values(30, 10.0); 

	        x[0] = raw_type(x1); 
	        x[1] = raw_type(x2); 
	        
	        VecCreateNest(PETSC_COMM_WORLD, 2, NULL, x, & raw_type(merged_vec)); 
	        DVectord result = merged_mat * merged_vec; 

	        bool verbose = false; 
	        EigenValueSlover<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> slepc; 

	        slepc.portion_of_spectrum("smallest_real"); 
	        slepc.number_of_eigenvalues(1); 
	        slepc.solver_type("arnoldi"); 

	        slepc.verbose(verbose); 
	        slepc.tol(1e-12); 
	        slepc.solve(merged_mat); 
	        slepc.print_eigenpairs(); 


	        VecDestroy(& raw_type(merged_vec)); 
	        MatDestroy(& raw_type(merged_mat)); 
	    }



	    void tr_more_sorensen_eigen_test()
	    {
	    	if(mpi_world_size() != 1)
	    		return; 

	       	DVectord x_w1  = values(4, 10);
			DVectord expected_woods = values(4, 1);
			Woods<DMatrixd, DVectord> fun_woods;

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
			subproblem->atol(1e-14); 
			subproblem->rtol(1e-14); 

			TrustRegion<DMatrixd, DVectord> tr_solver(subproblem);
			tr_solver.verbose(verbose);
			tr_solver.max_it(500); 
			tr_solver.atol(1e-11); 
			tr_solver.rtol(1e-14); 
			tr_solver.stol(1e-12); 
			tr_solver.solve(fun_woods, x_w1);				

			auto eigen_solver = std::make_shared<EigenValueSlover<DMatrixd, DVectord, PETSC_EXPERIMENTAL> >();
			eigen_solver->solver_type("arpack");
			
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
