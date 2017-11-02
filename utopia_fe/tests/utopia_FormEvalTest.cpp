#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"

namespace utopia {

	static void run_interp_vec_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		auto mesh = std::make_shared<Mesh>();
		mesh->make_triangle();

		auto Vx = FunctionSpaceT(mesh);
		auto Vy = FunctionSpaceT(mesh);
		auto V = Vx * Vy;

		auto u = trial(V);
		auto v = test(V);

		DVectord c_uk = values(6, 1.);
		{
			Write<DVectord> w(c_uk);
			c_uk.set(0, 1.0);
			c_uk.set(1, 2.0);
			c_uk.set(2, 3.0);
			c_uk.set(3, 1.0);
			c_uk.set(4, 2.0);
			c_uk.set(5, 3.0);
		}

		auto uk = interpolate(c_uk, u);
		auto b_form = integral( inner(grad(uk) * grad(u), grad(v)) );

		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(b_form);

		ElementMatrix mat;

		FormEvaluator<Backend> eval;
		eval.eval(b_form, mat, ctx, true);
		disp(mat);
	}

	static void run_interp_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		auto mesh = std::make_shared<Mesh>();
		mesh->make_triangle();

		auto V  = FunctionSpaceT(mesh);

		auto u = trial(V);
		auto v = test(V);

		DVectord c_uk = values(3, 1.);
		{
			Write<DVectord> w(c_uk);
			c_uk.set(0, 1.0);
			c_uk.set(1, 2.0);
			c_uk.set(2, 3.0);
		}

		auto uk = interpolate(c_uk, u);
		auto b_form = integral(inner(uk * grad(u), grad(v)) );

		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(b_form);

		ElementMatrix mat;

		FormEvaluator<Backend> eval;
		eval.eval(b_form, mat, ctx, true);
		disp(mat);
	}

	static void run_linear_elasticity()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		double mu = 1.;
		double lambda = 1.;

		// std::function<double(const AssemblyContext<Backend> &ctx)> fun_mu = [](const AssemblyContext<Backend> &ctx) -> double
		// {
		// 	return 1.;
		// };

		// auto mu = ctx_fun(fun_mu);

		auto mesh = std::make_shared<Mesh>();
		mesh->make_triangle();

		auto Vx = FunctionSpaceT(mesh);
		auto Vy = FunctionSpaceT(mesh);
		auto V  = Vx * Vy;

		auto u = trial(V);
		auto v = test(V);

		auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) ); 
		auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

		auto b_form = integral((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v)));

		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(b_form);

		ElementMatrix mat;

		FormEvaluator<Backend> eval;
		eval.eval(b_form, mat, ctx, true);
		disp(mat);
	}

	static void run_navier_stokes_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		auto mesh = std::make_shared<Mesh>();
		mesh->make_triangle();


		auto Vx = FunctionSpaceT(mesh);
		auto Vy = FunctionSpaceT(mesh);
		auto V  = Vx * Vy;

		auto Q =  FunctionSpaceT(mesh);

		auto u = trial(V);
		auto v = test(V);

		auto p = trial(Q);
		auto q = test(Q);

		const double mu = 1.;
		const double rho = 1.;
		
		const int n_dofs = 3 * 2;
		Vectord u_vector = values(n_dofs, 0.1);
		{
			Write<Vectord> w(u_vector);
			u_vector.set(1, 0.2);
		}

		auto uk = interpolate(u_vector, u);
		auto g_uk = grad(uk);

		auto e = mu * (transpose(grad(u)) + grad(u));
		auto b_form_11 = integral(inner(e, grad(v)) + rho * inner(g_uk * u, v));
		auto b_form_12 = integral(-inner(p, div(v)));
		auto b_form_21 = integral(inner(div(u), q));

		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(b_form_11);

		ElementMatrix mat;

		FormEvaluator<Backend> eval;
		eval.eval(b_form_11, mat, ctx, true);
		disp(mat);

		ctx.init_bilinear(b_form_12);
		eval.eval(b_form_12, mat, ctx, true);
		disp(mat);


		ctx.init_bilinear(b_form_21);
		eval.eval(b_form_21, mat, ctx, true);
		disp(mat);
	}

	static void run_vector_form_eval_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		auto mesh = std::make_shared<Mesh>();
		mesh->make_triangle();
		

		auto Vx = FunctionSpaceT(mesh);
		auto Vy = FunctionSpaceT(mesh);
		auto Vz = FunctionSpaceT(mesh);

		auto V = Vx * Vy * Vz;
		
		auto u = trial(V);
		auto v = test(V);

		Matrixd A = identity(3, 3);

		auto mass = integral( dot(u, v) );
		auto laplacian = integral( dot(transpose(A) * grad(u), grad(v)) );

		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(mass);

		ElementMatrix mat;

		FormEvaluator<Backend> eval;
		eval.eval(mass, mat, ctx, true);
		disp(mat);

		eval.eval(laplacian, mat, ctx, true);
		disp(mat);
	}

	static void run_scalar_form_eval_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		static const int Backend = Traits<FunctionSpaceT>::Backend;

		const int order = 2;

		auto mesh = std::make_shared<Mesh>();
		mesh->make_triangle(order);

		auto V = FunctionSpaceT(mesh, order);
		auto u = trial(V);
		auto v = test(V);

		Matrixd A = identity(2, 2);
		{
			Write<Matrixd> w(A);
		 	A.set(0, 0, 0.1);
		}

		auto mass         = integral(dot(u, v));
		auto diff_op    = 0.1 * (-abs(integral(1. * dot( (A  + transpose(A)) * grad(u), grad(v)) - 0.1 * mass, 0))) + 0.9 * integral( dot(2. * u, v), 2);
		auto linear_form  = integral(0.5 * dot(coeff(0.1), v));

		static_assert( (IsSubTree<TrialFunction<utopia::Any>, decltype(mass)>::value), "could not find function" );
		static_assert( (IsSubTree<TestFunction<utopia::Any>,  decltype(mass)>::value), "could not find function" );

		static_assert( (IsSubTree<TrialFunction<utopia::Any>, decltype(diff_op)>::value), "could not find function" );
		static_assert( (IsSubTree<TestFunction<utopia::Any>,  decltype(diff_op)>::value), "could not find function" );

		static_assert( (IsSubTree<TestFunction<utopia::Any>,   decltype(linear_form)>::value), "could not find function" );
		static_assert( !(IsSubTree<TrialFunction<utopia::Any>, decltype(linear_form)>::value), "should not find function" );

		ElementMatrix mat;
		ElementVector vec;

		FormEvaluator<Backend> eval;

		//init context only once for highest order assembly
		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(mass);

		eval.eval(diff_op, mat, ctx, true);

		disp("----------------------");
		disp(mat);
		disp("----------------------");

		eval.eval(mass, mat, ctx, true);
		disp(mat);
		disp("----------------------");

		ctx.init_linear(linear_form);
		eval.eval(linear_form, vec, ctx, true);
		disp(vec);
		disp("----------------------");
	}

	static void run_mixed_form_eval_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		auto mesh = std::make_shared<Mesh>();
		mesh->make_triangle();

		//scalar
		auto U = FunctionSpaceT(mesh);
		
		//vector
		auto Vx = FunctionSpaceT(mesh);
		auto Vy = FunctionSpaceT(mesh);
		auto V  = Vx * Vy;

		auto u = trial(U);
		auto v = test(V);

		Matrixd A = identity(2, 2);
		{
			Write<Matrixd> w(A);
		 	A.set(1, 1, 2.);
		}

		auto mixed = integral( dot(grad(u), 0.1 * A * v) );

		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(mixed);

		ElementMatrix mat;

		FormEvaluator<Backend> eval;
		eval.eval(mixed, mat, ctx, true);
		disp(mat);
	}

	void run_form_eval_test(libMesh::LibMeshInit &init)
	{
		std::cout << "run_scalar_form_eval_test:" << std::endl;
		run_scalar_form_eval_test();
		std::cout << "run_vector_form_eval_test:" << std::endl;
		run_vector_form_eval_test();
		std::cout << "run_mixed_form_eval_test:" << std::endl;
		run_mixed_form_eval_test();
		std::cout << "run_linear_elasticity:" << std::endl;
		run_linear_elasticity();
		std::cout << "run_interp_test:" << std::endl;
		run_interp_test();
		std::cout << "run_interp_vec_test:" << std::endl;
		run_interp_vec_test();
		std::cout << "run_navier_stokes_test:" << std::endl;
		run_navier_stokes_test();
	}
}
