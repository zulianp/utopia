#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"


namespace utopia {

	template<class MeshT, class FunctionSpaceT>
	class FormEvalTest {
	public:
		typedef utopia::Traits<FunctionSpaceT> TraitsT;
		static const int Backend = TraitsT::Backend;

		static void run_all_on(const std::shared_ptr<MeshT> &mesh)
		{
			std::cout << "--------------- [run_scalar_form_sum_eval_test] ---------------:" << std::endl;
			run_scalar_form_sum_eval_test(mesh);
			std::cout << "--------------- [run_scalar_form_eval_test] ---------------:" << std::endl;
			run_scalar_form_eval_test(mesh);
			std::cout << "--------------- [run_vector_form_eval_test] ---------------:" << std::endl;
			run_vector_form_eval_test(mesh);
			std::cout << "--------------- [run_mixed_form_eval_test] ---------------:" << std::endl;
			run_mixed_form_eval_test(mesh);
			std::cout << "--------------- [run_linear_elasticity] ---------------:" << std::endl;
			run_linear_elasticity(mesh);
			std::cout << "--------------- [run_interp_test] ---------------:" << std::endl;
			run_interp_test(mesh);
			std::cout << "--------------- [run_interp_vec_test] ---------------:" << std::endl;
			run_interp_vec_test(mesh);
			std::cout << "--------------- [run_navier_stokes_test] ---------------:" << std::endl;
			run_navier_stokes_test(mesh);
			std::cout << "--------------- [leastsquares_helmoholtz] ---------------:" << std::endl;
			leastsquares_helmoholtz(mesh);
			std::cout << "--------------- [run_eq_form_eval_test] ---------------:" << std::endl;
			run_eq_form_eval_test(mesh);
		}


		template<class Form>
		static void assemble_bilinear_and_print(const Form &form)
		{
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
			AssemblyContext<Backend> ctx;
			ctx.init_linear(form);

			ElementVector vec;

			FormEvaluator<Backend> eval;
			eval.eval(form, vec, ctx, true);

			std::cout << tree_format(form.getClass()) << std::endl;
			disp(vec);
		}

		static void run_interp_vec_test(const std::shared_ptr<MeshT> &mesh)
		{
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

			assemble_bilinear_and_print(b_form);
		}

		static void run_interp_test(const std::shared_ptr<MeshT> &mesh)
		{
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

			assemble_bilinear_and_print(b_form);
		}

		static void run_linear_elasticity(const std::shared_ptr<MeshT> &mesh)
		{
			double mu = 1.;
			double lambda = 1.;

			auto Vx = FunctionSpaceT(mesh);
			auto Vy = FunctionSpaceT(mesh);
			auto V  = Vx * Vy;

			auto u = trial(V);
			auto v = test(V);

			auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) ); 
			auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

			auto b_form = integral((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v)));

			assemble_bilinear_and_print(b_form);
		}

		static void run_navier_stokes_test(const std::shared_ptr<MeshT> &mesh)
		{
			auto Vx = FunctionSpaceT(mesh);
			auto Vy = FunctionSpaceT(mesh);
			auto V  = Vx * Vy;

			auto Q =  FunctionSpaceT(mesh);

			auto u = trial(V);
			auto v = test(V);

			auto p = trial(Q);
			auto q = test(Q);

			const double mu  = 1.;
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
			auto b_form_12 = -integral(inner(p, div(v)));
			auto b_form_21 = integral(inner(div(u), q));

			assemble_bilinear_and_print(b_form_11);
			assemble_bilinear_and_print(b_form_12);
			assemble_bilinear_and_print(b_form_21);
		}

		static void run_vector_form_eval_test(const std::shared_ptr<MeshT> &mesh)
		{
			auto Vx = FunctionSpaceT(mesh);
			auto Vy = FunctionSpaceT(mesh);
			auto Vz = FunctionSpaceT(mesh);

			auto V = Vx * Vy * Vz;

			auto u = trial(V);
			auto v = test(V);

			Matrixd A = identity(3, 3);

			auto mass = integral( inner(u, v) );
			auto laplacian = integral( inner(transpose(A) * grad(u), grad(v)) );

			assemble_bilinear_and_print(mass);
			assemble_bilinear_and_print(laplacian);
		}

		static void run_scalar_form_eval_test(const std::shared_ptr<MeshT> &mesh)
		{
			const int order = 1;
			auto V = FunctionSpaceT(mesh, order);
			auto u = trial(V);
			auto v = test(V);

			Matrixd A = identity(2, 2);
			{
				Write<Matrixd> w(A);
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


		static void run_scalar_form_sum_eval_test(const std::shared_ptr<MeshT> &mesh)
		{
			const int order = 1;

			auto V = FunctionSpaceT(mesh, order);
			auto u = trial(V);
			auto v = test(V);

			Matrixd A = identity(2, 2);
			{
				Write<Matrixd> w(A);
				A.set(0, 0, 0.1);
			}

			auto diff_op = integral(inner(u, v) + inner(u, v));

			assemble_bilinear_and_print(diff_op);
		}

		static void run_mixed_form_eval_test(const std::shared_ptr<MeshT> &mesh)
		{
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

			auto mixed = integral( inner(grad(u), 0.1 * A * v) );
			assemble_bilinear_and_print(mixed);
		}

		static void leastsquares_helmoholtz(const std::shared_ptr<MeshT> &mesh)
		{
			//scalar
			auto U = FunctionSpaceT(mesh);

			//vector
			auto Vx = FunctionSpaceT(mesh);
			auto Vy = FunctionSpaceT(mesh);
			auto V  = Vx * Vy;

			auto u = trial(U);
			auto v = test(U);	

			auto s = trial(V);
			auto q = test(V);

			// strong_enforce( boundary_conditions(u == coeff(0.0), {0, 1, 2, 3}) );

			//bilinear forms
			double c = -100.0;
			double beta = 0.99;

			auto eq_11 = integral((c*c) * inner(u, v) + inner(grad(u), grad(v)));
			auto eq_12 = integral(c * inner(div(s), v) + inner(s, grad(v)));
			auto eq_21 = integral(c * inner(u, div(q)) + inner(grad(u), q));
			auto eq_22 = integral(inner(s, q) + inner(div(s), div(q)) + beta * inner(curl(s), curl(q)));


			//linear forms
			auto f = coeff(1);
			auto rhs_1 = integral(c * dot(f, v));
			auto rhs_2 = integral(dot(f, div(u)));

			assemble_bilinear_and_print(eq_11);
			assemble_bilinear_and_print(eq_12);
			assemble_bilinear_and_print(eq_21);
			assemble_bilinear_and_print(eq_22);	
		}

		static void run_eq_form_eval_test(const std::shared_ptr<MeshT> &mesh)
		{
			auto V = FunctionSpaceT(mesh);

			auto u = trial(V);
			auto v = test(V);
			auto bilinear_form = integral(inner(grad(u), grad(v)));
			auto linear_form   = integral(inner(coeff(1.), v));
			auto eq = bilinear_form == linear_form;

			ElementMatrix mat;
			ElementVector vec;
			FormEvaluator<Backend> eval;
			eval.eval_equation(eq, mat, vec, true);

			disp(mat);
			disp(vec);
		}

		static void run_time_form_eval_test(const std::shared_ptr<MeshT> &mesh)
		{
			auto V = FunctionSpaceT(mesh);

			auto u = trial(V);
			auto v = test(V);

			DVectord c_uk = zeros(3);
			auto uk = interpolate(c_uk, u);
			auto eq = integral(inner(dt(uk), v)) == integral(inner(coeff(1.), v) -  0.1 * inner(grad(uk), grad(v)));
		}
	};

	void run_form_eval_test(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<utopia::Mesh>();
		mesh->make_triangle();

		FormEvalTest<utopia::Mesh, utopia::HMFESpace>  test;
		test.run_all_on(mesh);
	}
}
