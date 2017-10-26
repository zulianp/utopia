#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"

namespace utopia {

	void run_vector_form_eval_test()
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
		auto laplacian = integral( dot(grad(u), grad(v)) );

		AssemblyContext<Backend> ctx;
		ctx.init_bilinear(mass);

		ElementMatrix mat;

		FormEvaluator<Backend> eval;
		eval.eval(mass, mat, ctx, true);
		disp(mat);

		eval.eval(laplacian, mat, ctx, true);
		disp(mat);
	}

	void run_scalar_form_eval_test()
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

	void run_mixed_form_eval_test()
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

		auto mixed = integral( dot(grad(u), v) );

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
	}
}
