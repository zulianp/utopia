#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"

namespace utopia {

	static void run_linear_elasticity()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		double mu = 1.;
		double lambda = 1.;

		auto Vx = FunctionSpaceT();
		auto Vy = FunctionSpaceT();
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


		auto Vx = FunctionSpaceT();
		auto Vy = FunctionSpaceT();
		auto V  = Vx * Vy;

		auto Q =  FunctionSpaceT();

		auto u = trial(V);
		auto v = test(V);

		auto p = trial(Q);
		auto q = test(Q);

		const double mu = 1.;
		const double rho = 1.;
		auto inc_cond = integral(inner(q, div(u)));
		
		// const int n_dofs = 3 * 2;
		// DVector u_vector = zeros(n_dofs);
		// auto uk = interpolate(u, u_vector);
		// auto g_uk = grad(uk);
		// auto e = 0.5 * (transpose(grad(u))+ grad(u));
		// auto mom1 =  integral(inner(2. * mu * e, grad(v)) + rho * inner(g_uk * u, v);
		// auto mom2 = -integral(inner(p, div(v)));
 			
	}

	static void run_vector_form_eval_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;
		

		auto Vx = FunctionSpaceT();
		auto Vy = FunctionSpaceT();
		auto Vz = FunctionSpaceT();

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

		auto V = FunctionSpaceT();
		auto u = trial(V);
		auto v = test(V);

		Matrixd A = identity(2, 2);
		{
			Write<Matrixd> w(A);
		 	A.set(0, 0, 0.1);
		}

		auto mass         = integral(dot(u, v));
		auto laplacian    = 0.1 * (-abs(integral(1. * dot( (A  + transpose(A)) * grad(u), grad(v)) - 0.1 * mass, 0))) + 0.9 * integral( dot(2. * u, v), 2);
		auto linear_form  = integral(0.5 * dot(coeff(0.1), v));

		ElementMatrix mat;
		ElementVector vec;

		FormEvaluator<Backend> eval;

		//init context only once for highest order assembly
		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(mass);

		eval.eval(laplacian, mat, ctx, true);

		disp(mat);

		eval.eval(mass, mat, ctx, true);
		disp(mat);

		ctx.init_linear(linear_form);
		eval.eval(linear_form, vec, ctx, true);
		disp(vec);
	}

	static void run_mixed_form_eval_test()
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		typedef utopia::Traits<HMFESpace> TraitsT;
		static const int Backend = TraitsT::Backend;

		//scalar
		auto U = FunctionSpaceT();
		
		//vector
		auto Vx = FunctionSpaceT();
		auto Vy = FunctionSpaceT();
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
		run_scalar_form_eval_test();
		run_vector_form_eval_test();
		run_mixed_form_eval_test();
		run_linear_elasticity();
	}
}
