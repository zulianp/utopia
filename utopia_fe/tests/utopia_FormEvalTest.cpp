#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"

namespace utopia {
	void run_form_eval_test(libMesh::LibMeshInit &init)
	{
		//scalar function spaces
		auto V1 = FunctionSpace<EmptyFunctionSpace>();
		auto V2 = FunctionSpace<EmptyFunctionSpace>();
		auto V3 = FunctionSpace<EmptyFunctionSpace>();

		//vector function space
		auto V_vec = V1 * V2;
		auto V = kron_prod(V_vec, V3);
		auto &V_vec_ref    = V.get<0>();
		auto &V_scalar_ref = V.get<1>();

		//vector fe
		auto uv = trial(V_vec_ref);
		auto vv = test(V_vec_ref);

		//scalar fe
		auto us = trial(V_scalar_ref);
		auto vs = test(V_scalar_ref);		

		Vectord r = values(2, 1.);

		//forms
		auto bilinear_form = integral( dot(grad(us), vv) );
		auto linear_form   = integral( dot(r, vs) );

		FormEvaluator<decltype(bilinear_form), EMPTY_BACKEND_FLAG> eval_bilinear;
		FormEvaluator<decltype(linear_form),   EMPTY_BACKEND_FLAG> eval_linear;

		Matrixd mat;
		eval_bilinear.eval(bilinear_form, mat, true);

		Vectord vec;
		eval_linear.eval(linear_form, vec, true);
	}
}
