#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"

#include "utopia_LibMeshBackend.hpp"
#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"

#include "libmesh/exodusII_io.h"
#include <algorithm>

namespace utopia {

	/*
		TODO list
		- boundary conditions
		- temporal derivative
		- variational inequality
		- projection
		- blocks/side-sets/node-sets functions
		- local-2-global
		- multi-linear form
	*/


	class EqPrinter {
	public:
		template<class Eq>
		void operator()(const int index, const Eq &eq) const {
			std::cout << "equation: " << index << std::endl;
			std::cout << tree_format(eq.getClass()) << std::endl;
		}
	};

	


	template<class SpaceInput, class FunctionSpaceT>
	class FormEvalTest {
	public:
		typedef utopia::Traits<FunctionSpaceT> TraitsT;
		typedef typename TraitsT::Matrix ElementMatrix;
		typedef typename TraitsT::Vector ElementVector;

		static const int Backend = TraitsT::Backend;

		static void run_all_on(const std::shared_ptr<SpaceInput> &space_input)
		{
			run_linear_form_test(space_input);
			run_bilinear_form_test(space_input);
			run_scalar_form_sum_eval_test(space_input);
			run_scalar_form_eval_test(space_input);
			run_vector_form_eval_test(space_input);
			run_mixed_form_eval_test(space_input);
			run_linear_elasticity(space_input);
			run_interp_test(space_input);
			run_interp_vec_test(space_input);
			run_navier_stokes_test(space_input);
			leastsquares_helmoholtz(space_input);
			run_eq_form_eval_test(space_input);
			run_local_2_global_test(space_input);
		}

		template<class Equations>
		static void init_context(Equations &eqs)
		{

		}

		template<class Form>
		static void assemble_bilinear_and_print(const Form &form)
		{
			static_assert(IsForm<Form>::value,      "must be a form");
			static_assert(IsForm<Form>::order == 2, "must be a bilinear form");
			static_assert(IsForm<Form>::has_trial,  "must have trial function");
			static_assert(IsForm<Form>::has_test,   "must have test function");
			static_assert(IsForm<Form>::has_fun,    "must have function");

			AssemblyContext<Backend> ctx;
			ctx.init_bilinear(form);

			ElementMatrix mat;

			FormEvaluator<Backend> eval;
			eval.eval(form, mat, ctx, true);

			std::cout << tree_format(form.getClass()) << std::endl;
			disp(mat);
		}


		template<class Form>
		static void assemble_linear_and_print(const Form &form)
		{
			static_assert(IsForm<Form>::value,      "must be a form");
			static_assert(IsForm<Form>::order == 1, "must be a linear form");
			static_assert(IsForm<Form>::has_test,   "must be a linear form");
			static_assert(IsForm<Form>::has_fun,    "must have function");

			AssemblyContext<Backend> ctx;
			ctx.init_linear(form);

			ElementVector vec;

			FormEvaluator<Backend> eval;
			eval.eval(form, vec, ctx, true);

			std::cout << tree_format(form.getClass()) << std::endl;
			disp(vec);
		}

		static void run_interp_vec_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_interp_vec_test]" << std::endl;

			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
			auto V = Vx * Vy;

			auto u = trial(V);
			auto v = test(V);

			ElementVector c_uk = values(27, 1.);
			{
				Write<ElementVector> w(c_uk);
				c_uk.set(0, 1.0);
				c_uk.set(1, 2.0);
				c_uk.set(2, 3.0);
				c_uk.set(3, 1.0);
				c_uk.set(4, 2.0);
				c_uk.set(5, 3.0);
			}

			auto uk = interpolate(c_uk, u);
			auto b_form = integral( inner(grad(uk) * grad(u), grad(v)) );

			assemble_bilinear_and_print(b_form);
		}

		static void run_interp_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_interp_test]" << std::endl;

			auto V  = FunctionSpaceT(space_input);

			auto u = trial(V);
			auto v = test(V);

			V.initialize();


			//large enough vector
			ElementVector c_uk = values(27, 1.);
			{
				Write<ElementVector> w(c_uk);
				c_uk.set(0, 1.0);
				c_uk.set(1, 2.0);
				c_uk.set(2, 3.0);
			}

			auto uk = interpolate(c_uk, u);
			auto b_form = integral(inner(uk * grad(u), grad(v)) );

			assemble_bilinear_and_print(b_form);
		}

		static void run_linear_elasticity(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_linear_elasticity]" << std::endl;

			double mu = 1.;
			double lambda = 1.;

			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
			auto V  = Vx * Vy;

			auto u = trial(V);
			auto v = test(V);

			auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) ); 
			auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

			auto b_form = integral((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v)));

			assemble_bilinear_and_print(b_form);
		}


		static void run_nonlinear_elasticity(const std::shared_ptr<SpaceInput> &space_input)
		{

			std::cout << "[run_nonlinear_elasticity]" << std::endl;

			double mu = 1.;
			double lambda = 1.;

			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
			auto V  = Vx * Vy;

			auto u = trial(V);
			auto v = test(V);

			DVectord sol = zeros(27);
			each_write(sol, [](const SizeType i) -> double {
				return 0.1 + i/27.;
			});

			auto uk = interpolate(sol, u);

			auto F 		 = identity() + grad(uk);
			auto F_inv   = inv(F);
			auto F_inv_t = transpose(F_inv);
			auto g_uk    = grad(uk);

			auto J = det(F);
			auto C = F + transpose(F);
			auto t = trace(C);
			
			auto P = mu * (F - F_inv_t) + lambda * logn(J) * F_inv_t;
			//compressible neo-hookean
			auto l_form = inner(P, grad(v)) * dX;
			// auto b_form = (
			// 				mu * inner(grad(u), grad(v))
			// 	 		    - (lambda * logn(J) - mu) * inner(transpose(F_inv * grad(u)), F_inv * grad(v))
			// 	 		    + lambda * inner(F_inv_t, grad(u)) * inner(F_inv_t, grad(v))
			// 	 		  ) * dX;

			auto b_form = (inner(P, grad(u)) * inner(P, grad(v))) * dX;
		


			static_assert(IsForm<decltype(l_form)>::value,      "must be a form");
			static_assert(IsForm<decltype(l_form)>::order == 1, "must be a form");


			assemble_linear_and_print(l_form);
			assemble_bilinear_and_print(b_form);
		}

		static void run_navier_stokes_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_navier_stokes_test]" << std::endl;

			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
			auto V  = Vx * Vy;

			auto Q =  FunctionSpaceT(space_input);

			auto u = trial(V);
			auto v = test(V);

			auto p = trial(Q);
			auto q = test(Q);

			const double mu  = 1.;
			const double rho = 1.;

			//2 elements with 3 nodes and 2+1 dofs per node
			const int n_dofs = 2 * 3 * (2 + 1);
			ElementVector u_vector = values(n_dofs, 0.1);
			{
				Write<ElementVector> w(u_vector);
				u_vector.set(1, 0.2);
			}

			auto uk = interpolate(u_vector, u);
			auto g_uk = grad(uk);

			auto e = mu * (transpose(grad(u)) + grad(u));
			auto b_form_11 = integral(inner(e, grad(v)) + rho * inner(g_uk * u, v));
			auto b_form_12 = -integral(inner(p, div(v)));
			auto b_form_21 = integral(inner(div(u), q));


			ElementVector r_v = values(2, 1.);
			auto l_form_1 =  integral(inner(coeff(1.), q));
			auto l_form_2 =  integral(inner(coeff(r_v), v));

			auto lhs = b_form_11 + b_form_12 + b_form_21;
			auto rhs = l_form_1 + l_form_2;

			assemble_bilinear_and_print(lhs);
			assemble_linear_and_print(rhs);


			static_assert(IsSubTree<Integral<utopia::Any>, decltype(lhs)>::value, "must be an integral");
			static_assert(IsSubTree<Interpolate<Any, Any>, decltype(lhs)>::value, "non-linear term not detected");
		}

		static void run_vector_form_eval_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_vector_form_eval_test]" << std::endl;

			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
			auto Vz = FunctionSpaceT(space_input);

			auto V = Vx * Vy * Vz;

			auto u = trial(V);
			auto v = test(V);

			ElementMatrix A = identity(3, 3);

			auto mass = integral( inner(u, v) );
			auto laplacian = integral( inner(transpose(A) * grad(u), grad(v)) );

			assemble_bilinear_and_print(mass);
			assemble_bilinear_and_print(laplacian);
		}

		static void run_scalar_form_eval_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_scalar_form_eval_test]" << std::endl;

			
			auto V = FunctionSpaceT(space_input);
			auto u = trial(V);
			auto v = test(V);

			ElementMatrix A = identity(2, 2);
			{
				Write<ElementMatrix> w(A);
				A.set(0, 0, 0.1);
			}

			auto mass        = integral(inner(u, v));
			auto diff_op     = 0.1 * (-abs(integral(1. * inner( (A  + transpose(A)) * grad(u), grad(v)) , 0) - 0.1 * mass)) + 0.9 * integral( inner(2. * u, v), 2);
			auto linear_form = integral(0.5 * inner(coeff(0.1), v));

			static_assert( (IsSubTree<TrialFunction<utopia::Any>,  decltype(mass)>::value), 	"could not find function" );
			static_assert( (IsSubTree<TestFunction<utopia::Any>,   decltype(mass)>::value), 	"could not find function" );

			static_assert( (IsSubTree<TrialFunction<utopia::Any>,  decltype(diff_op)>::value), 	"could not find function" );
			static_assert( (IsSubTree<TestFunction<utopia::Any>,   decltype(diff_op)>::value), 	"could not find function" );

			static_assert( (IsSubTree<TestFunction<utopia::Any>,   decltype(linear_form)>::value), "could not find function" );
			static_assert( !(IsSubTree<TrialFunction<utopia::Any>, decltype(linear_form)>::value), "should not find function" );

			assemble_bilinear_and_print(mass);
			assemble_bilinear_and_print(diff_op);
			assemble_linear_and_print(linear_form);
		}


		static void run_scalar_form_sum_eval_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_scalar_form_sum_eval_test]" << std::endl;

			auto V = FunctionSpaceT(space_input);
			auto u = trial(V);
			auto v = test(V);

			ElementMatrix A = identity(2, 2);
			{
				Write<ElementMatrix> w(A);
				A.set(0, 0, 0.1);
			}

			auto diff_op = integral(inner(u, v) + inner(u, v));



			assemble_bilinear_and_print(diff_op);
		}

		static void run_mixed_form_eval_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_mixed_form_eval_test]" << std::endl;

			//scalar
			auto U = FunctionSpaceT(space_input);

			//vector
			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
			auto V  = Vx * Vy;

			auto u = trial(U);
			auto v = test(V);
			auto w = test(U);

			ElementMatrix A = identity(2, 2);
			{
				Write<ElementMatrix> w(A);
				A.set(1, 1, 2.);
			}

			auto mixed = integral( inner(grad(u), 0.1 * A * v) ) + integral(inner(grad(u), grad(w)));
			assemble_bilinear_and_print(mixed);
		}

		static void leastsquares_helmoholtz(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[leastsquares_helmoholtz]" << std::endl;

			//scalar
			auto V = FunctionSpaceT(space_input);

			//vector
			auto Qx = FunctionSpaceT(space_input);
			auto Qy = FunctionSpaceT(space_input);
			auto Q  = Qx * Qy;

			auto u = trial(V);
			auto v = test(V);	

			auto s = trial(Q);
			auto q = test(Q);

			// strong_enforce( boundary_conditions(u == coeff(0.0), {0, 1, 2, 3}) );

			//bilinear forms
			double c = -100.0;
			double beta = 0.99;

			auto eq_11 = integral((c*c) * inner(u, v) + inner(grad(u), grad(v)));
			auto eq_12 = integral(c * inner(div(s), v) + inner(s, grad(v)));
			auto eq_21 = integral(c * inner(u, div(q)) + inner(grad(u), q));
			auto eq_22 = integral(inner(s, q) + inner(div(s), div(q)) + beta * inner(curl(s), curl(q)));


			//linear forms
			auto rhs_1 = integral(c * inner(coeff(1.), v));
			auto rhs_2 = integral(inner(coeff(1.), div(q)));

			assemble_bilinear_and_print(eq_11 + eq_12 + eq_21 + eq_22);
			assemble_linear_and_print(rhs_1 + rhs_2);
		}

		static void run_eq_form_eval_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_eq_form_eval_test]" << std::endl;

			auto V = FunctionSpaceT(space_input);

			auto u = trial(V);
			auto v = test(V);
			auto bilinear_form = integral(inner(grad(u), grad(v)));
			auto linear_form   = integral(inner(coeff(1.), v));
			auto eq = bilinear_form == linear_form;

			ElementMatrix mat;
			ElementVector vec;
			FormEvaluator<Backend> eval;
			eval.eval_equation(eq, mat, vec, true);

			// disp(mat);
			// disp(vec);
		}

		static void run_time_form_eval_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_time_form_eval_test]" << std::endl;

			auto V = FunctionSpaceT(space_input);

			auto u = trial(V);
			auto v = test(V);

			ElementVector c_uk = zeros(3);
			auto uk = interpolate(c_uk, u);
			auto eq = integral(inner(dt(uk), v)) == integral(inner(coeff(1.), v) -  0.1 * inner(grad(uk), grad(v)));
		}

		static void run_mixed_fe_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_mixed_fe_test]" << std::endl;

			auto V  = FunctionSpaceT(space_input);
			auto Wx = FunctionSpaceT(space_input);
			auto Wy = FunctionSpaceT(space_input);
			auto W  = Wx * Wy;

			auto Q = FunctionSpaceT(space_input);
			// auto T = TensorFunctionSpace(space_input);
			// auto T = VectorFunctionSpace(space_input);

			//create a mixed function space
			auto M = mixed(V, W, Q);

			//create trial functions
			auto trial_funs = trials(M);
			auto test_funs  = tests(M);

			//get individual functions
			auto u = std::get<0>(trial_funs);
			auto z = std::get<1>(trial_funs);

			auto v = std::get<0>(test_funs);
			auto w = std::get<1>(test_funs);
		}

		static void run_linear_form_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_linear_form_test]" << std::endl;

			auto V = FunctionSpaceT(space_input);
			auto u = trial(V);
			auto v = test(V);

			auto linear_form = integral(0.5 * inner(coeff(0.1), v));
			assemble_linear_and_print(linear_form);
		}

		static void run_bilinear_form_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_bilinear_form_test]" << std::endl;


			auto V = FunctionSpaceT(space_input);
			auto u = trial(V);
			auto v = test(V);

			auto bilinear_form = integral(inner(u, v));
			assemble_bilinear_and_print(bilinear_form);
		}



		static void run_nonlinear_form_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_nonlinear_form_test]" << std::endl;


			auto V = FunctionSpaceT(space_input);
			auto u = trial(V);
			auto v = test(V);


			DVectord sol;
			auto uk = interpolate(sol, u);
			auto alpha = coeff(0.5);
			auto R = u * alpha * (coeff(1.) - uk) * (uk - alpha);
			// auto R = u * alpha * (coeff(3.) - uk);
			// auto R = u * (uk - alpha);
			// auto R = u * uk;

			V.initialize();
			sol = local_values(V.dof_map().n_local_dofs(), 2.);

			auto bilinear_form = integral(inner(R, v));
			// assemble_bilinear_and_print(bilinear_form);

			auto linear_form = integral(inner(uk * alpha * (coeff(1.) - uk) * (uk - alpha), v));
			assemble_linear_and_print(linear_form);
		}

		static void run_local_2_global_test(const std::shared_ptr<SpaceInput> &space_input)
		{

			// std::cout << "[run_local_2_global_test]" << std::endl;

			// auto V = FunctionSpaceT(space_input);

			// auto u = trial(V);
			// auto v = test(V);
			// auto bilinear_form = integral(inner(grad(u), grad(v)));
			// auto linear_form   = integral(inner(coeff(1.), v));
			// auto eq = bilinear_form == linear_form;

			// DSMatrixd mat;
			// DVectord vec;
			// assemble_expression_v<FunctionSpaceT>(eq, mat, vec);
		}
	};

	typedef FormEvalTest<libMesh::EquationSystems, utopia::LibMeshFunctionSpace> LibMeshFormEvalTest;

	template<typename F>
	static void run_libmesh_test(libMesh::LibMeshInit &init, F fun)
	{	
		auto lm_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		
		const unsigned int n = 10;
		libMesh::MeshTools::Generation::build_square(*lm_mesh,
			n, n,
			0, 1,
			0, 1.,
			libMesh::QUAD8);

		std::shared_ptr<libMesh::EquationSystems> equation_systems;
		LibMeshFormEvalTest lm_test;

		equation_systems = std::make_shared<libMesh::EquationSystems>(*lm_mesh);
		fun(lm_test, equation_systems);
	}

	void run_libmesh_eval_test(libMesh::LibMeshInit &init)
	{
		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_linear_form_test");
		// 	lm_test.run_linear_form_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_bilinear_form_test");
		// 	lm_test.run_bilinear_form_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_scalar_form_sum_eval_test");
		// 	lm_test.run_scalar_form_sum_eval_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_scalar_form_eval_test");
		// 	lm_test.run_scalar_form_eval_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_vector_form_eval_test");
		// 	lm_test.run_vector_form_eval_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_mixed_form_eval_test");
		// 	lm_test.run_mixed_form_eval_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_linear_elasticity");
		// 	lm_test.run_linear_elasticity(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_interp_test");
		// 	lm_test.run_interp_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_interp_vec_test");
		// 	lm_test.run_interp_vec_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_navier_stokes_test");
		// 	lm_test.run_navier_stokes_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("leastsquares_helmoholtz");
		// 	lm_test.leastsquares_helmoholtz(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_eq_form_eval_test");
		// 	lm_test.run_eq_form_eval_test(es);
		// });

		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_local_2_global_test");
		// 	lm_test.run_local_2_global_test(es);
		// });


		// run_libmesh_test(init,[](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_nonlinear_form_test");
		// 	lm_test.run_nonlinear_form_test(es);
		// });

		run_libmesh_test(init,[](
			LibMeshFormEvalTest &lm_test,
			const std::shared_ptr<libMesh::EquationSystems> &es) {
			es->add_system<libMesh::LinearImplicitSystem>("run_nonlinear_elasticity");
			lm_test.run_nonlinear_elasticity(es);
		});

		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	auto &sys = es->add_system<libMesh::LinearImplicitSystem>("test_equations");

		// 	//space for u
		// 	auto V = LibMeshFunctionSpace(es);

		// 	auto u = trial(V);
		// 	auto v = test(V);
		// 	DVectord sol;
		// 	const bool success = solve(
		// 		equations(
		// 			inner(grad(u), grad(v)) * dX == inner(coeff(0.), v) * dX
		// 		),
		// 		constraints(
		// 			boundary_conditions(u == coeff(-0.2), {1}),
		// 			boundary_conditions(u == coeff(0.2),  {2})
		// 		),
		// 		sol);


		// 	//go back to libmesh
		// 	convert(sol, *sys.solution);
		// 	sys.solution->close();
		// 	libMesh::ExodusII_IO(V.mesh()).write_equation_systems ("test_equations.e", *es);
		// });


		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	auto &sys = es->add_system<libMesh::LinearImplicitSystem>("test_elasticity");

		// 	//space for u
		// 	auto Vx = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::FIRST, "disp_x");
		// 	auto Vy = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::FIRST, "disp_y");
		// 	auto V = Vx * Vy;

		// 	auto u = trial(V);
		// 	auto v = test(V);

		// 	auto ux = u[0];
		// 	auto uy = u[1];

		// 	const double mu = 1;
		// 	const double lambda = 1;

		// 	auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) ); 
		// 	auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

		// 	LMDenseVector z = zeros(2);

		// 	DVectord sol;
		// 	const bool success = solve(
		// 		equations(
		// 			((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v))) * dX == inner(coeff(z), v) * dX

		// 		),
		// 		constraints(
		// 			boundary_conditions(uy == coeff(0.2),  {0}),
		// 			boundary_conditions(uy == coeff(-0.2), {2}),
		// 			boundary_conditions(ux == coeff(0.0),  {0, 2})
		// 		),
		// 		sol);


		// 	//go back to libmesh
		// 	convert(sol, *sys.solution);
		// 	sys.solution->close();
		// 	libMesh::ExodusII_IO(Vx.mesh()).write_equation_systems ("test_elasticity.e", *es);
		// });


		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	auto &sys = es->add_system<libMesh::LinearImplicitSystem>("non_linear_laplacian");

		// 	//space for u
		// 	auto V = LibMeshFunctionSpace(es);

		// 	auto u = trial(V);
		// 	auto v = test(V);
			
		// 	DVectord sol;
		// 	auto uk = interpolate(sol, u);

		// 	if(nl_solve(
		// 		equations(
		// 			(inner(grad(u), grad(v)) + inner(grad(uk) * u, grad(v))) * dX == inner(coeff(0.0), v) * dX
		// 		),
		// 		constraints(
		// 			boundary_conditions(u == coeff(-0.2), {1}),
		// 			boundary_conditions(u == coeff(0.2),  {2})
		// 		),
		// 		sol)) {
		// 		//go back to libmesh
		// 		convert(sol, *sys.solution);
		// 		sys.solution->close();
		// 		libMesh::ExodusII_IO(V.mesh()).write_equation_systems ("non_linear_laplacian.e", *es);
		// 	} else {
		// 		std::cerr << "[Error] solver failed to converge" << std::endl;
		// 	}
		// });


		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	auto &sys = es->add_system<libMesh::LinearImplicitSystem>("navier_stokes");

	
		// 	auto Vx = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::SECOND, "vel_x");
		// 	auto Vy = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::SECOND, "vel_y");
		// 	auto V  = Vx * Vy;

		// 	auto Q =  LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::FIRST, "pressure");

		// 	auto u = trial(V);
		// 	auto v = test(V);

		// 	auto ux = u[0];
		// 	auto uy = u[1];

		// 	auto p = trial(Q);
		// 	auto q = test(Q);




		// 	const double mu  = 1.;
		// 	const double rho = 1.;
		// 	const double dt = 0.01;
		// 	const std::size_t n_ts = 10;

		// 	DVectord sol;
		// 	DVectord sol_old;


		// 	auto uk = interpolate(sol, u);
		// 	auto uk_old = interpolate(sol_old, u);

		// 	auto g_uk = grad(uk);

		// 	auto e         = mu * (transpose(grad(u)) + grad(u));
		// 	auto b_form_11 = integral(inner(e, grad(v)) + rho * inner(g_uk * u, v));
		// 	auto b_form_12 = integral(inner(p, div(v)));
		// 	auto b_form_21 = integral(inner(div(u), q));

		// 	LMDenseVector r_v = zeros(2);
		// 	auto l_form_1 =  integral(inner(coeff(0.), q));
		// 	auto l_form_2 =  integral(inner(uk, v));

		// 	auto lhs = b_form_11 + b_form_12 + b_form_21;
		// 	auto rhs = l_form_1  + l_form_2;

		// 	if(nl_solve(
		// 		// nl_implicit_euler(
		// 		equations(
		// 			lhs == rhs
		// 		),
		// 		constraints(
		// 			boundary_conditions(ux == coeff(0.), {0, 1, 3}),
		// 			boundary_conditions(ux == coeff(0.1), {2}),
		// 			boundary_conditions(uy == coeff(0.), {0, 1, 2, 3})
		// 		),
		// 		// sol_old,
		// 		sol
		// 		//, n_ts
		// 		)) {
		// 		//go back to libmesh
		// 		convert(sol, *sys.solution);
		// 		sys.solution->close();
		// 		libMesh::ExodusII_IO(Q.mesh()).write_equation_systems ("navier_stokes.e", *es);
		// 	} else {
		// 		std::cerr << "[Error] solver failed to converge" << std::endl;
		// 	}
		// });



		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	auto &sys = es->add_system<libMesh::LinearImplicitSystem>("reaction_diffusion");
		// 	auto V = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::FIRST, "u");
		// 	auto u = trial(V);
		// 	auto v = test(V);

		// 	const double dt = 0.01;
		// 	const std::size_t n_ts = 40;

		// 	DVectord sol;
		// 	DVectord sol_old;

		// 	//if_else(cond, val_if, val_else)

		// 	// auto uk     = interpolate(sol, in_block(u == 0.1, {1, 2}) || 
		// 								   // in_block(u == 0.,  {0, 3}) );

		// 	auto uk     = interpolate(sol, u);
		// 	auto uk_old = interpolate(sol_old, u);

		// 	auto a = 100.;
		// 	auto alpha = coeff(0.5);
		// 	auto R = uk * a * (coeff(1.) - uk) * (uk - alpha);

		// 	auto f_rhs = ctx_fun< std::vector<double> >([&u](const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
		// 		const auto &pts = ctx.fe()[0]->get_xyz();

		// 		const auto n = pts.size();
		// 		std::vector<double> ret(n);
				
		// 		for(std::size_t i = 0; i != n; ++i) {
		// 			double x = (pts[i](0) - 0.5);
		// 			double y = (pts[i](1) - 0.5);
		// 			double dist = std::sqrt(x*x + y*y);
		// 			if(dist < 0.2)
		// 				ret[i] = 2;
		// 			else
		// 				ret[i] = 0.;
		// 		}

		// 	 	return ret;
		// 	});

		// 	if(nl_implicit_euler(
		// 		equations(
		// 			( inner(u, v) + inner( dt * grad(u), grad(v)) ) * dX == ( dt * inner(R, v) + inner(uk_old, v) + dt * inner(f_rhs, v)) * dX
		// 		),
		// 		constraints(
		// 			boundary_conditions(u == coeff(0.),  {1, 3}),
		// 			boundary_conditions(u == coeff(0.),  {0}),
		// 			boundary_conditions(u == coeff(0.0), {2})
		// 		),
		// 		sol_old,
		// 		sol, 
		// 		dt,
		// 		n_ts
		// 		)) 
		// 	{
		// 		//post process
		// 	} else {
		// 		std::cerr << "[Error] solver failed to converge" << std::endl;
		// 	}
		// });


		/////////////////////////////////////////////////////////////////////

		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	auto &sys = es->add_system<libMesh::LinearImplicitSystem>("test_user_kernel");

		// 	//space for u
		// 	auto V = LibMeshFunctionSpace(es);

		// 	auto u = trial(V);
		// 	auto v = test(V);
		// 	DVectord sol;
		// 	const bool success = solve(
		// 		equations(
		// 			bilinear_kernel<double>(
		// 				[](libMesh::FEBase &trial,
		// 				   libMesh::FEBase &test, 
		// 				   const unsigned int trial_index,
		// 				   const unsigned int test_index,
		// 				   const unsigned int qp) -> libMesh::Real {

		// 				return (trial.get_dphi()[trial_index][qp] * test.get_dphi()[test_index][qp]) * trial.get_JxW()[qp];
					
		// 			}, 0) == inner(coeff(0.), v) * dX
		// 		),
		// 		constraints(
		// 			boundary_conditions(u == coeff(-0.2), {1}),
		// 			boundary_conditions(u == coeff(0.2),  {2})
		// 		),
		// 		sol);


		// 	//go back to libmesh
		// 	convert(sol, *sys.solution);
		// 	sys.solution->close();
		// 	libMesh::ExodusII_IO(V.mesh()).write_equation_systems ("test_user_kernel.e", *es);
		// });


		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	auto &sys = es->add_system<libMesh::LinearImplicitSystem>("transient_navier_stokes");

		
		// 	auto Vx = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::SECOND, "vel_x");
		// 	auto Vy = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::SECOND, "vel_y");
		// 	auto V  = Vx * Vy;

		// 	auto Q =  LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::FIRST, "pressure");

		// 	auto u = trial(V);
		// 	auto v = test(V);

		// 	auto ux = u[0];
		// 	auto uy = u[1];

		// 	auto p = trial(Q);
		// 	auto q = test(Q);


		// 	const double mu  = 1.;
		// 	const double rho = 1.;
		// 	const double dt = 0.005;
		// 	const std::size_t n_ts = 40;

		// 	DVectord sol;
		// 	DVectord sol_old;


		// 	auto uk 	= interpolate(sol, u);
		// 	auto uk_old = interpolate(sol_old, u);

		// 	auto g_uk = grad(uk);

		// 	auto e         = mu * (transpose(grad(u)) + grad(u));
		// 	auto b_form_11 = integral(inner(u, v)) + dt * integral(inner(e, grad(v)) - (rho * dt) * inner(g_uk * u, v));
		// 	auto b_form_12 = dt * integral(inner(p, div(v)));
		// 	auto b_form_21 = integral(inner(div(u), q));

		// 	LMDenseVector r_v = zeros(2);
		// 	auto l_form_1 =  integral(inner(coeff(0.), q));
		// 	auto l_form_2 =  integral(inner(uk_old, v)); //zero forcing term

		// 	auto lhs = b_form_11 + b_form_12 + b_form_21;
		// 	auto rhs = l_form_1  + l_form_2;

		// 	if(nl_implicit_euler(
		// 		equations(
		// 			lhs == rhs
		// 		),
		// 		constraints(
		// 			// boundary_conditions(ux == coeff(0.), {0, 1, 3}),
		// 			// boundary_conditions(ux == coeff(0.1), {2}),
		// 			// boundary_conditions(uy == coeff(0.), {0, 1, 2, 3})
		// 			boundary_conditions(uy == coeff(0.),   {0, 1, 2, 3}),
		// 			boundary_conditions(ux == coeff(0.),   {0, 2}),
		// 			boundary_conditions(ux == coeff(0.1),  {1, 3})
		// 		),
		// 		sol_old,
		// 		sol,
		// 		dt,
		// 		n_ts
		// 		)) {
		// 		//....
		// 	} else {
		// 		std::cerr << "[Error] solver failed to converge" << std::endl;
		// 	}
		// });

	}

	void run_form_eval_test(libMesh::LibMeshInit &init)
	{

		/////////////////////////////////////////////////////////////////////
		//HOMEMADE
		// auto mesh = std::make_shared<utopia::Mesh>();
		// mesh->make_triangle();
		// FormEvalTest<utopia::Mesh, utopia::HMFESpace>  test;
		// test.run_all_on(mesh);

		/////////////////////////////////////////////////////////////////////


		/////////////////////////////////////////////////////////////////////
		//LIBMESH
		run_libmesh_eval_test(init);

	}
}
