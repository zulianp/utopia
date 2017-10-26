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

		auto laplacian    = integral( dot(grad(u), grad(v)) );
		auto mass         = integral( dot(u, v) );
		auto linear_form  = integral( dot(coeff(0.1), v) );

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
