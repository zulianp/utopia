#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"

namespace utopia {
	void run_form_eval_test(libMesh::LibMeshInit &init)
	{
		typedef utopia::HMFESpace FunctionSpaceT;
		static const int Backend = Traits<FunctionSpaceT>::Backend;

		auto V = FunctionSpaceT();
		auto u = trial(V);
		auto v = test(V);

		auto laplacian = integral( dot(grad(u), grad(v)) );
		auto mass = integral( dot(u, v) );
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

	// void run_mixed_form_eval_test()
	// {
	// 	typedef utopia::HMFESpace FunctionSpaceT;
	// 	// typedef utopia::EmptyFunctionSpace FunctionSpaceT;

	// 	//scalar function spaces
	// 	auto V1 = FunctionSpaceT();
	// 	auto V2 = FunctionSpaceT();
	// 	auto V3 = FunctionSpaceT();

	// 	//vector function space
	// 	auto V_vec = V1 * V2;
	// 	auto V = kron_prod(V_vec, V3);
	// 	auto &V_vec_ref    = V.get<0>();
	// 	auto &V_scalar_ref = V.get<1>();

	// 	//vector fe
	// 	auto uv = trial(V_vec_ref);
	// 	auto vv = test(V_vec_ref);

	// 	//scalar fe
	// 	auto us = trial(V_scalar_ref);
	// 	auto vs = test(V_scalar_ref);		

	// 	Vectord r = values(2, 1.);

	// 	//forms
	// 	auto bilinear_form = integral( dot(grad(us), vv) );
	// 	auto linear_form   = integral( dot(r, vs) );

	// 	Matrixd mat;
	// 	Vectord vec;

	// 	FormEvaluator<HOMEMADE> eval;
	// 	eval.eval_bilinear(bilinear_form, mat, true);
	// 	eval.eval_linear(linear_form, vec, true);
	// }
}
