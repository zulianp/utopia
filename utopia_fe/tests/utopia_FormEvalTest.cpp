#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_function.hpp"
#include "utopia.hpp"

namespace utopia {
	void run_form_eval_test(libMesh::LibMeshInit &init)
	{
		//scalar function spaces
		auto V1 = FunctionSpace<EmptyFunctionSpace>();
		auto V2 = FunctionSpace<EmptyFunctionSpace>();

		//vector function space
		auto V = V1 * V2;

		auto u = trial(V);
		auto v = test(V);

		Vectord r = values(2, 1.);

		auto bilinear_form = integral(dot(grad(u), grad(v)));
		auto linear_form   = integral(dot(r, v));

		FormEvaluator<decltype(bilinear_form), EMPTY_BACKEND_FLAG> eval_bilinear;

		Matrixd mat;
		eval_bilinear.eval(bilinear_form, mat, true);


		FormEvaluator<decltype(linear_form), EMPTY_BACKEND_FLAG> eval_linear;

		Vectord vec;
		eval_linear.eval(linear_form, vec, true);
	}
}