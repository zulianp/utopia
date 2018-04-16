/*
* @Author: kopanicakova
* @Date:   2018-02-06 17:47:26
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-04-10 11:00:58
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
	class PetscNonlinearSolverTest {
	public:
		
		void run()
		{
			UTOPIA_RUN_TEST(petsc_ngs_test);
			UTOPIA_RUN_TEST(petsc_gss_newton_test);
			UTOPIA_RUN_TEST(petsc_newton_test);
			UTOPIA_RUN_TEST(petsc_newton_rosenbrock_test);
			UTOPIA_RUN_TEST(petsc_sparse_semismooth_newton_test);
			UTOPIA_RUN_TEST(petsc_sparse_nonlinear_semismooth_newton_test);
			UTOPIA_RUN_TEST(petsc_direct_solver_newton_test);
			UTOPIA_RUN_TEST(petsc_newton_test_out_info);
			UTOPIA_RUN_TEST(petsc_sparse_newton_test);
			UTOPIA_RUN_TEST(petsc_newton_petsc_cg_test);
			UTOPIA_RUN_TEST(petsc_tr_rr_test);
			UTOPIA_RUN_TEST(petsc_mprgp_test);
			UTOPIA_RUN_TEST(petsc_inexact_newton_test);
			UTOPIA_RUN_TEST(petsc_snes_test); 
			UTOPIA_RUN_TEST(petsc_sparse_newton_snes_test); 
		}

		void petsc_ngs_test()
		{
			const int n = 50;
			DSMatrixd m = sparse(n, n, 3);
			assemble_laplacian_1D(n, m);
			// const double ub = 100.0;
			const double ub = 1.;
			const bool use_line_search = true;
			// const bool use_line_search = false;
			const int max_it = n * 60;
			const bool verbose = false;
			const int n_local_sweeps = 3;

			{
				Range r = row_range(m);
				Write<DSMatrixd> w(m);
				if(r.begin() == 0) {
					m.set(0, 0, 1.);
					m.set(0, 1, 0);
				}

				if(r.end() == n) {
					m.set(n-1, n-1, 1.);
					m.set(n-1, n-2, 0);
				}
			}

			DVectord rhs = values(n, 1.);
			rhs *= 1./(n-1);

			{ 
			    //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
				Range r = range(rhs);
				Write<DVectord> w(rhs);

				if(r.begin() == 0) {
					rhs.set(0, 0);
				}

				if(r.end() == n) {
					rhs.set(n-1, 0.);
				}
			}

			DVectord upper_bound = values(n, ub);
			DVectord solution    = zeros(n);

			ProjectedGaussSeidel<DSMatrixd, DVectord, -1> pgs;
			pgs.max_it(max_it);
			pgs.verbose(verbose);
			pgs.set_use_line_search(use_line_search);
			pgs.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

			Chrono c;
			c.start();
			
			pgs.solve(m, rhs, solution);
			
			c.stop();
			// if(mpi_world_rank() == 0) std::cout << c << std::endl;


			DVectord solution_u = zeros(n);
			ProjectedGaussSeidel<DSMatrixd, DVectord, -1> pgs_u;
			pgs_u.verbose(verbose);
			pgs_u.set_n_local_sweeps(n_local_sweeps);
			pgs_u.set_use_line_search(use_line_search);
			pgs_u.set_use_symmetric_sweep(true);
			pgs_u.max_it(max_it);

			pgs_u.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

			c.start();
			
			pgs_u.solve(m, rhs, solution_u);
			
			c.stop();
			// if(mpi_world_rank() == 0) std::cout << c << std::endl;

			double diff = norm2(solution_u - solution);
			double res_norm = norm2(m * solution_u - rhs);

			// disp(res_norm);

			if(diff > 1e-5) {
				std::cerr << "[Error] different implementations of pgs gives different results, diff: " << diff << std::endl;
			}

			assert(approxeq(solution_u, solution, 1e-5));

			//standard gs with MatSOR
			// GaussSeidel<DSMatrixd, DVectord> gs;
			// gs.verbose(verbose);
			// gs.max_it(max_it);
			// // gs.sweeps(1);
			// solution.set(0.);
			
			// c.start();
			
			// gs.solve(m, rhs, solution);

			// c.stop();

			// if(mpi_world_rank() == 0) std::cout << c << std::endl;

			// double res_norm_ref = norm2(m * solution - rhs);

			// // disp(res_norm_ref);

			// if(diff > 1e-5) {
			// 	std::cerr << "[Error] different implementations of pgs gives different results, diff: " << diff << std::endl;
			// }

			// assert(approxeq(solution_u, solution, 1e-5));
		}

		void petsc_gss_newton_test()
		{
			typedef std::function<void(const DSMatrixd &, const DVectord &, const DVectord &, DVectord &, DVectord &)> F;

			const int n = mpi_world_size() * 4;
			DVectord sol  = zeros(n);
			DVectord upbo = values(n, 1.);
			DSMatrixd A   = identity(n, n);
			DVectord rhs  = values(n, 3.);

			DVectord lambda, d;
			F f = [&lambda, &d, &upbo](const DSMatrixd &H, const DVectord &g, const DVectord &x, DVectord &active, DVectord &value) 
			{
				lambda = (upbo - H * x);
				d = lambda + (x - upbo);	

				Read<DVectord> r_d(d);
				Read<DVectord> r_u(upbo);
				Write<DVectord> w_d(active);
				Write<DVectord> w_v(value);

				auto rr = range(x);
				for (SizeType i = rr.begin(); i != rr.end(); i++) {
					if (d.get(i) >= -1e-16) {
						active.set(i, 1.0);
						value.set(i, upbo.get(i));
					} else {
						active.set(i, 0.0);
						value.set(i,  0.);
					}
				}
			};

			auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
			GenericSemismoothNewton<DSMatrixd, DVectord, F> solver(f, linear_solver);

			solver.solve(A, rhs, sol);
		}
		
		
		void petsc_mprgp_test()
		{
			const SizeType n = 50;
			const PetscScalar h = 1.0/(n-1);
			
			
			DSMatrixd A = sparse(n, n, 3);
			DVectord b, u, l;
			
			// 1d laplace
			{
				Write<DSMatrixd> w_A(A);
				const Range r = row_range(A);
				
				for(SizeType i = r.begin(); i != r.end(); ++i) {
					if(i > 0) {
						A.add(i, i - 1, -1.0);
					}
					
					if(i < n-1) {
						A.add(i, i + 1, -1.0);
					}
					
					A.add(i, i, 2.0);
				}
			}
			
			A = (n-1)*A;
			
			// bc conditions
			{
				Range r = row_range(A);
				Write<DSMatrixd> w(A);
				if(r.begin() == 0) {
					A.set(0, 0, 1.);
					A.set(0, 1, 0);
				}
				
				if(r.end() == n) {
					A.set(n-1, n-1, 1.);
					A.set(n-1, n-2, 0.);
				}
			}
			
			l = values(n, -1.);
			u = values(n, 1.);
			
			b = values(n, 50.);
			
			{
				Write<DVectord> w (b);
				Range rhs_range = range(b);;
				for (SizeType i = rhs_range.begin(); i != rhs_range.end() ; i++) {
					if(i > n/2) {
						b.set(i, -50);
					}
					
					if(i ==0 || i == n-1) {
						b.set(i, 0);
					}
				}
			}
			
			b *= h;
			
			// just testing
			// u =  values(local_size(u).get(0), 999);
			// l =  values(local_size(u).get(0), -999);
			
			auto box = make_box_constaints(make_ref(l), make_ref(u));
			
			//FIXME the following unilateral constraints do not work with MPRGP
			// auto box = make_upper_bound_constraints(make_ref(u));
			
			MPRGP<DSMatrixd, DVectord> mprgp;
			mprgp.set_box_constraints(make_ref(box));
			
			// initial guess
			DVectord x = 0 * b;
			
			mprgp.verbose(false);
			mprgp.solve(A, b, x);
			
			
			DVectord x_0 = 0. * x;

			auto lsolver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>();
			// auto lsolver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
			// auto lsolver = std::make_shared<GMRES<DSMatrixd, DVectord>>();
			// auto lsolver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord, HOMEMADE>>();
			lsolver->atol(1e-15);
			lsolver->rtol(1e-15);
			lsolver->stol(1e-15);
			// lsolver->verbose(true);

			SemismoothNewton<DSMatrixd, DVectord> nlsolver(lsolver);
			nlsolver.set_box_constraints(box);
			// nlsolver.verbose(true);

			nlsolver.max_it(200);
			nlsolver.solve(A, b, x_0);		

			// disp(l);
			// disp(u);	
			assert(approxeq(x, x_0));
		}
		
		void petsc_tr_rr_test()
		{
			// rosenbrock test
			if(mpi_world_size() == 1)
			{
				DVectord x = values(10, 2);
				TestFunctionND_1<DMatrixd, DVectord> fun2(x.size().get(0));
				DVectord expected = values(x.size().get(0), 0.468919);
				
				Parameters params;
				params.atol(1e-10);
				params.rtol(1e-10);
				params.stol(1e-10);
				params.verbose(false);
				
				auto subproblem = std::make_shared<utopia::KSP_TR<DMatrixd, DVectord> >();
				TrustRegion<DMatrixd, DVectord> tr_solver(subproblem);
				tr_solver.set_parameters(params);
				tr_solver.solve(fun2, x);
				
				params.trust_region_alg(STEIHAUG_TOINT_TAG);
				x = values(10, 2);
				trust_region_solve(fun2, x, params);
				
				params.trust_region_alg(NASH_TAG);
				x = values(10, 2);
				trust_region_solve(fun2, x, params);
				
				params.trust_region_alg(LANCZOS_TAG);
				x = values(10, 2);
				trust_region_solve(fun2, x, params);
				
				params.trust_region_alg(CGNE_TAG);
				x = values(10, 2);
				trust_region_solve(fun2, x, params);
				
			}
		}

		
		void petsc_newton_test_out_info()
		{
			if(mpi_world_size() > 10) return;
			
			auto lsolver = std::make_shared< BiCGStab<DMatrixd, DVectord> >();
			Newton<DMatrixd, DVectord> nlsolver(lsolver);
			
			Parameters params;
			params.verbose(false);
			params.linear_solver_verbose(false);
			nlsolver.set_parameters(params);
			
			DVectord x = values(10, 2);
			TestFunctionND_1<DMatrixd, DVectord> fun2(x.size().get(0));
			
			DVectord expected = values(x.size().get(0), 0.468919);
			nlsolver.solve(fun2, x);
			assert(approxeq(expected, x));
		}
		
		void petsc_sparse_newton_test()
		{
			if(mpi_world_size() > 10) return;
			
			auto lsolver = std::make_shared< BiCGStab<DSMatrixd, DVectord> >();
			Newton<DSMatrixd, DVectord> nlsolver(lsolver);
			nlsolver.enable_differentiation_control(false);
			
			Parameters params;
			params.verbose(false);
			params.linear_solver_verbose(false);
			nlsolver.set_parameters(params);
			
			SimpleQuadraticFunction<DSMatrixd, DVectord> fun;
			
			DVectord x = values(10, 2.);
			DVectord expected = zeros(x.size());
			
			nlsolver.solve(fun, x);
			assert(approxeq(expected, x));
		}
		
		void petsc_newton_test()
		{
			if(mpi_world_size() > 10) return;
			
			auto lsolver = std::make_shared< BiCGStab<DMatrixd, DVectord> >();
			Newton<DMatrixd, DVectord> nlsolver(lsolver);
			nlsolver.enable_differentiation_control(false);
			
			Parameters params;
			params.atol(1e-15);
			params.rtol(1e-15);
			params.stol(1e-15);
			params.verbose(false);
			nlsolver.set_parameters(params);
			
			SimpleQuadraticFunction<DMatrixd, DVectord> fun;
			
			DVectord x = values(10, 2.);
			DVectord expected = zeros(x.size());
			
			nlsolver.solve(fun, x);
			assert(approxeq(expected, x));
			
			x = values(10, 2.0);
			TestFunctionND_1<DMatrixd, DVectord> fun2(x.size().get(0));
			
			expected = values(x.size().get(0), 0.468919);
			nlsolver.solve(fun2, x);
			assert(approxeq(expected, x));
		}

		void petsc_inexact_newton_test()
		{
			if(mpi_world_size() > 10) return;
			
			Parameters params;
			params.atol(1e-15);
			params.rtol(1e-15);
			params.stol(1e-15);
			params.verbose(false);
			
			auto lsolver = std::make_shared< BiCGStab<DMatrixd, DVectord> >();
			InexactNewton<DMatrixd, DVectord> nlsolver(lsolver);
			nlsolver.set_parameters(params);
			
			auto hess_approx_BFGS   = std::make_shared<BFGS<DMatrixd, DVectord> >();
			nlsolver.set_hessian_approximation_strategy(hess_approx_BFGS);
			
			
			SimpleQuadraticFunction<DMatrixd, DVectord> fun;
			
			DVectord x = values(10, 2.);
			DVectord expected_1 = zeros(x.size());
			
			nlsolver.solve(fun, x);
			assert(approxeq(expected_1, x));
			
			TestFunctionND_1<DMatrixd, DVectord> fun2(x.size().get(0));
			x = values(10, 2.0);
			DVectord expected_2 = values(x.size().get(0), 0.468919);
			nlsolver.solve(fun2, x);
			assert(approxeq(expected_2, x));
			
			// -------------------------------------- SR1 test ------------------
			auto hess_approx_SR1    = std::make_shared<SR1<DMatrixd, DVectord> >();
			nlsolver.set_hessian_approximation_strategy(hess_approx_SR1);
			
			x = values(10, 2.);
			nlsolver.solve(fun, x);
			assert(approxeq(expected_1, x));
			
			x = values(10, 2.0);
			nlsolver.solve(fun2, x);
			assert(approxeq(expected_2, x));			
		}

		void petsc_newton_rosenbrock_test()
		{
			auto lsolver = std::make_shared< BiCGStab<DMatrixd, DVectord> >();
			Newton<DMatrixd, DVectord> nlsolver(lsolver);
			nlsolver.enable_differentiation_control(false);
			
			Parameters params;
			params.rtol(1e-15);
			params.stol(1e-15);
			params.atol(1e-15);
			params.verbose(false);
			nlsolver.set_parameters(params);
			
			DVectord expected_rosenbrock;
			DVectord x0;
			
			if (mpi_world_size() <= 2) {
				RosenbrockGeneric<DMatrixd, DVectord> r_generic_2d;
				expected_rosenbrock = values(2, 1.0);
				x0 = values(2, 2.0);
				nlsolver.solve(r_generic_2d, x0);
				assert(approxeq(expected_rosenbrock, x0));
			}
			
			if (mpi_world_size() <= 3) {
				RosenbrockGeneric<DMatrixd, DVectord> r_generic_3d;
				expected_rosenbrock = values(3, 1.0);
				x0 = values(3, -2.0);
				nlsolver.solve(r_generic_3d, x0);
				assert(approxeq(expected_rosenbrock, x0));
			}
			
			if (mpi_world_size() <= 6) {
				RosenbrockGeneric<DMatrixd, DVectord> r_generic_6d;
				expected_rosenbrock = values(6, 1.0);
				x0 = values(6, 2.0);
				nlsolver.solve(r_generic_6d, x0);
				assert(approxeq(expected_rosenbrock, x0));
			}
		}
		
		void petsc_sparse_semismooth_newton_test()
		{
			auto lsolver = std::make_shared<Factorization<DSMatrixd, DVectord>>();

			DSMatrixd A;
			DVectord b, ub;
			
			SemismoothNewton<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> petsc_ss_newton(lsolver);
			SemismoothNewton<DSMatrixd, DVectord, HOMEMADE> homemade_ss_newton(lsolver);
			
			// initial guess
			DVectord x_0 = values(_n, 0.0);
			
			ExampleTestCase2<DSMatrixd, DVectord> example;
			example.getOperators(_n, A, b, ub);

			const double scale_factor = 1;
			A *= scale_factor;
			b *= scale_factor;
			ub *= scale_factor;
			
			auto box = make_upper_bound_constraints(make_ref(ub));
			petsc_ss_newton.set_box_constraints(box);
			petsc_ss_newton.atol(1e-18);
			petsc_ss_newton.rtol(1e-15);
			petsc_ss_newton.stol(1e-16);
			petsc_ss_newton.max_it(400);
			petsc_ss_newton.solve(A, b, x_0);
		

			DVectord hm_x_0 = values(_n, 0.0);
			homemade_ss_newton.set_box_constraints(box);
			homemade_ss_newton.stol(1e-16);
			// homemade_ss_newton.verbose(true);
			homemade_ss_newton.solve(A, b, hm_x_0);


			// x_0.implementation().set_name("x");
			// hm_x_0.implementation().set_name("y");
		

			x_0 *= 1./scale_factor;
			hm_x_0 *= 1./scale_factor;

			// write("x_p.m", x_0);
			// write("x_u.m", hm_x_0);

			if(!approxeq(x_0, hm_x_0, 1e-14)) {
				DVectord diff = hm_x_0 - x_0;
			// 	disp(diff);
				double norm_diff = norm2(diff)/double(norm2(hm_x_0));
				std::cout << "norm_diff: " << norm_diff << std::endl;
			}

			assert(approxeq(x_0, hm_x_0, 1e-8));
		}
		
		void petsc_sparse_nonlinear_semismooth_newton_test()
		{
			auto lsolver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>();
			
			NonlinSemismoothNewton<DSMatrixd, DVectord> nlsolver(lsolver);
			nlsolver.enable_differentiation_control(false);
			Parameters params;
			params.verbose(false);
			nlsolver.set_parameters(params);
			
			DSMatrixd A, B;
			DVectord upbo;
			
			ExampleTestCase<DSMatrixd, DVectord> example;
			example.getOperators(_n, A, B, upbo);
			
			DVectord rhs = values(_n, 60);
			{
				Write<DVectord> w(rhs);
				Range rhs_range = range(rhs);
				if(rhs_range.begin() == 0) rhs.set(0, 0);
				if(rhs_range.end() == _n) rhs.set(_n - 1, 0);
			}
			
			QuadraticFunctionConstrained<DSMatrixd, DVectord> funn(rhs, A, B, upbo);
			
			auto box = make_upper_bound_constraints(make_ref(upbo));
			nlsolver.set_box_constraints(make_ref(box));
			
			nlsolver.solve(funn, rhs);
		}
		
		void petsc_direct_solver_newton_test()
		{
			auto lsolver = std::make_shared< Factorization<DSMatrixd, DVectord> >();
			
#ifdef PETSC_HAVE_MUMPS
			lsolver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
#endif //PETSC_HAVE_MUMPS
			
			Newton<DSMatrixd, DVectord> nlsolver(lsolver);
			nlsolver.enable_differentiation_control(false);
			
			Parameters params;
			params.verbose(false);
			params.linear_solver_verbose(false);
			nlsolver.set_parameters(params);
			
			SimpleQuadraticFunction<DSMatrixd, DVectord> fun;
			
			DVectord x = values(_n, 2.);
			DVectord expected = zeros(x.size());
			
			nlsolver.solve(fun, x);
			assert(approxeq(expected, x));
			
			
			auto lCG = std::make_shared< ConjugateGradient<DSMatrixd, DVectord> >();
			nlsolver.set_linear_solver(lCG);
			x = values(_n, 2.);
			nlsolver.solve(fun, x);
			assert(approxeq(expected, x));
		}
		
		
		void petsc_newton_petsc_cg_test()
		{
			using namespace std;
			
			
			//CG with diagonal preconditioner
			auto linear_solver  = make_shared< KSPSolver<DMatrixd, DVectord> >();
			auto preconditioner = make_shared< InvDiagPreconditioner<DMatrixd, DVectord> >();
			linear_solver->set_preconditioner(preconditioner);
			
			//Newton solver with cg linear solver
			Newton<DMatrixd, DVectord> newton_solver(linear_solver);
			newton_solver.verbose(false);
			
			const int n = 10;
			DVectord actual   = values(n, 2.);
			DVectord expected = values(n, 0.468919);
			
			TestFunctionND_1<DMatrixd, DVectord> fun(n);
			
			newton_solver.solve(fun, actual);
			assert(approxeq(expected, actual));
		}
		
		void petsc_newton_inexact_newton_with_KSP_test()
		{
			using namespace std;
			
			const bool verbose = false;
			
			auto linear_solver  = make_shared< KSPSolver<DMatrixd, DVectord> >();
			auto preconditioner = make_shared< InvDiagPreconditioner<DMatrixd, DVectord> >();
			linear_solver->set_preconditioner(preconditioner);
			
			
			InexactNewton<DMatrixd, DVectord> newton_solver(linear_solver);
			newton_solver.verbose(verbose);
			
			const int n = 10;
			DVectord actual   = values(n, 2.);
			DVectord expected = values(n, 0.468919);
			
			TestFunctionND_1<DMatrixd, DVectord> fun(n);
			
			newton_solver.solve(fun, actual);
			assert(approxeq(expected, actual));
		}

		void petsc_snes_test()
		{
			using namespace utopia;
			using namespace std;
		
			const static bool verbose = false;

			if(mpi_world_size() >= 10) return;
						
			auto linear_solver = make_shared< ConjugateGradient<DMatrixd, DVectord> >();

			SNESSolver<DMatrixd, DVectord, PETSC> nonlinear_solver(linear_solver); 
			nonlinear_solver.verbose(verbose); 

			if(mpi_world_size() == 1)
			{
				Rosenbrock<DMatrixd, DVectord> rosenbrock;
				DVectord expected_rosenbrock = values(2, 1.0);
				DVectord x0_ros   			 = values(2, 1.5);

				nonlinear_solver.solve(rosenbrock, x0_ros);


				expected_rosenbrock -= x0_ros; 
				double diff_rb = norm2(expected_rosenbrock);
				assert(approxeq(diff_rb, 0., 1e-6));


				// std::cout<<"--------------------------------------------------- \n"; 
				auto cg_home = std::make_shared<ConjugateGradient<DMatrixd, DVectord, HOMEMADE>>();
				cg_home->verbose(verbose); 

				SNESSolver<DMatrixd, DVectord, PETSC> nonlinear_solver2(cg_home); 
				nonlinear_solver2.verbose(verbose); 

				// reset IG  
				x0_ros   		    = values(2, 1.5);
				expected_rosenbrock = values(2, 1.0);
				nonlinear_solver2.solve(rosenbrock, x0_ros);


				expected_rosenbrock -= x0_ros; 
				diff_rb = norm2(expected_rosenbrock);
				assert(approxeq(diff_rb, 0., 1e-6));


				// std::cout<<"------------------ utopia-precond test --------------------------------- \n"; 

				auto preconditioner = make_shared< InvDiagPreconditioner<DMatrixd, DVectord> >();
				cg_home->set_preconditioner(preconditioner);

				SNESSolver<DMatrixd, DVectord, PETSC> nonlinear_solver3(cg_home); 
				nonlinear_solver3.verbose(verbose); 

			}
		}



		void petsc_sparse_newton_snes_test()
		{
			if(mpi_world_size() > 1) return;
			
			auto lsolver = std::make_shared< BiCGStab<DSMatrixd, DVectord> >();
			Newton<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> nlsolver(lsolver);
			
			Parameters params;
			params.verbose(false);
			params.linear_solver_verbose(false);
			nlsolver.set_parameters(params);
			
			SimpleQuadraticFunction<DSMatrixd, DVectord> fun;
			
			DVectord x = values(10, 2.);
			DVectord expected = zeros(x.size());
			
			nlsolver.solve(fun, x);
			assert(approxeq(expected, x));

			if(mpi_world_size() == 1)
			{
				Rosenbrock<DSMatrixd, DVectord> rosenbrock;
				DVectord expected_rosenbrock = values(2, 1.0);
				DVectord x0_ros   			 = values(2, 1.5);

				nlsolver.set_line_search_type("cp"); 
				nlsolver.set_line_search_order(3); 
				nlsolver.max_it(1000); 
				nlsolver.solve(rosenbrock, x0_ros); 

				expected_rosenbrock -= x0_ros; 
				double diff_rb = norm2(expected_rosenbrock);
				assert(approxeq(diff_rb, 0., 1e-6));
			}
		}


		PetscNonlinearSolverTest()
		: _n(100) { }
		
	private:
		int _n;
	};

#endif //WITH_PETSC
	

	void runPetscNonlinearSolversTest()
	{
		UTOPIA_UNIT_TEST_BEGIN("runPetscNonlinearSolverTest");
		#ifdef WITH_PETSC
				PetscNonlinearSolverTest().run();
		#endif
		UTOPIA_UNIT_TEST_END("runPetscNonlinearSolverTest");					
	}
}
