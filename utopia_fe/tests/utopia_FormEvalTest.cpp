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
#include "utopia_NonLinearFEFunction.hpp"
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

	//libmesh
	template<class FunctionSpaceT, class Expr, class GlobalMatrix, class GlobalVector>
	void element_assemble_expression_v(const libMesh::MeshBase::const_element_iterator &it, const Expr &expr, GlobalMatrix &mat, GlobalVector &vec)
	{
		typedef utopia::Traits<FunctionSpaceT> TraitsT;
		typedef typename TraitsT::Matrix ElementMatrix;
		typedef typename TraitsT::Vector ElementVector;

		static const int Backend = TraitsT::Backend;

		const auto &space = find_space<FunctionSpaceT>(expr);
		const auto &dof_map = space.dof_map();

		AssemblyContext<Backend> ctx;
		ctx.set_current_element((*it)->id());

		ElementMatrix el_mat;
		ElementVector el_vec;

		ctx.init_bilinear(expr);

		FormEvaluator<Backend> eval;
		eval.eval(expr, el_mat, el_vec, ctx, true);

		std::vector<libMesh::dof_id_type> dof_indices;
		dof_map.dof_indices(*it, dof_indices);

		dof_map.heterogenously_constrain_element_matrix_and_vector(el_mat.implementation(), el_vec.implementation(), dof_indices);

		add_matrix(el_mat.implementation(), dof_indices, dof_indices, mat);
		add_vector(el_vec.implementation(), dof_indices, vec);
	}

	//homemade
	template<class FunctionSpaceT, class Left, class Right, class Matrix, class Vector>
	void element_assemble_expression_v(const int element_index, const Equality<Left, Right> &equation, Matrix &mat, Vector &rhs)
	{
		std::cout << element_index << std::endl;
	}

	template<class FunctionSpaceT, class Left, class Right, class Matrix, class Vector>
	void assemble_expression_v(const Equality<Left, Right> &equation, Matrix &mat, Vector &vec)
	{
		auto &space = find_space<FunctionSpaceT>(equation);
		space.initialize();

		auto &m = space.mesh();
		auto &dof_map = space.dof_map();

		auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
								  *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

		mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
		vec = local_zeros(dof_map.n_local_dofs());

		Write<Matrix> w_m(mat);
		Write<Vector> w_v(vec);

		for(auto it = elements_begin(m); it != elements_end(m); ++it) {
			element_assemble_expression_v<FunctionSpaceT>(it, equation, mat, vec);
		}
	}

	template<class Matrix, class Vector, class Eqs>
	class NonLinearFEFunction : public Function<Matrix, Vector> {
	public:
		DEF_UTOPIA_SCALAR(Matrix)

		NonLinearFEFunction(const Eqs &eqs)
		: eqs_(eqs), first_(true)
		{}

		virtual ~NonLinearFEFunction() { }

		virtual bool value(const Vector &/*point*/, Scalar &/*value*/) const 
		{
			assert(false && "not implemented");
		    return 1.;
		}

		virtual bool gradient(const Vector &x, Vector &result) const
		{
			result = buff_mat * x - buff_vec;
		    return true;
		}


		virtual bool hessian(const Vector &, Matrix &H) const
		{
			H = buff_mat;
		    return true;
		}

		virtual bool update(const Vector &x) { 
			typedef decltype(eqs_.template get<0>()) Eq1;
			typedef typename FindFunctionSpace<Eq1>::Type FunctionSpaceT;
			auto &space = find_space<FunctionSpaceT>(eqs_);

			if(first_) {
				auto &dof_map = space.dof_map();
				auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
										  *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));
				
				buff_mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
				buff_vec = local_zeros(dof_map.n_local_dofs());
			} else {
				buff_mat *= 0.;
				buff_vec *= 0.;
			}

			{
				Write<DSMatrixd> w_m(buff_mat);
				Write<DVectord>  w_v(buff_vec);

				auto &m = space.mesh();

				for(auto it = elements_begin(m); it != elements_end(m); ++it) {
					element_assemble_expression_v<FunctionSpaceT>(it, eqs_, buff_mat, buff_vec);
				}
			}	

			return true;
		};

		bool first_;
		Eqs eqs_;

		Matrix buff_mat;
		Vector buff_vec;
	};


	template<class... Eqs, class... Constr>
	bool nl_solve(const Equations<Eqs...> &eqs, const FEConstraints<Constr...> &constr, DVectord &sol)
	{
		typedef typename GetFirst<Eqs...>::Type Eq1Type;
		typedef typename FindFunctionSpace<Eq1Type>::Type FunctionSpaceT;
		
		
		FEBackend<LIBMESH_TAG>::init_constraints(constr);

		auto &space = find_space<FunctionSpaceT>(eqs.template get<0>());
		space.initialize();

		sol = local_zeros(space.dof_map().n_local_dofs());

		NonLinearFEFunction<DSMatrixd, DVectord, Equations<Eqs...>> nl_fun(eqs);
		Newton<DSMatrixd, DVectord> solver(std::make_shared<Factorization<DSMatrixd, DVectord>>());
		solver.verbose(true);
		return solver.solve(nl_fun, sol);
	}


	template<class... Eqs, class... Constr>
	bool solve(const Equations<Eqs...> &eqs, const FEConstraints<Constr...> &constr, DVectord &sol)
	{

		//FIXME this stuff only works for the libmesh backend
		typedef typename GetFirst<Eqs...>::Type Eq1Type;
		typedef typename FindFunctionSpace<Eq1Type>::Type FunctionSpaceT;
		auto &space = find_space<FunctionSpaceT>(eqs.template get<0>());

		FEBackend<LIBMESH_TAG>::init_constraints(constr);


		
		space.initialize();
		auto &m = space.mesh();
		auto &dof_map = space.dof_map();
		auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
								  *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

		DSMatrixd mat;
		DVectord vec;

		mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
		vec = local_zeros(dof_map.n_local_dofs());

		{
			Write<DSMatrixd> w_m(mat);
			Write<DVectord>  w_v(vec);


			for(auto it = elements_begin(m); it != elements_end(m); ++it) {
				element_assemble_expression_v<FunctionSpaceT>(it, eqs, mat, vec);
			}
		}	


		// disp(mat);
		// disp(vec);

		sol = local_zeros(local_size(vec));
		Factorization<DSMatrixd, DVectord> solver;
		if(!solver.solve(mat, vec, sol)) {
			return false;
		}

		//if non-linear
		

		DVectord residual = mat * sol - vec;
		double norm_r = norm2(residual);
		std::cout << "norm_r: " << norm_r << std::endl;
		// disp(sol);

		// write("u_matrix_new.m", mat);
		return false;
	}



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
			AssemblyContext<Backend> ctx;
			ctx.init_bilinear(form);

			ElementMatrix mat;

			FormEvaluator<Backend> eval;
			eval.eval(form, mat, ctx, true);

			// std::cout << tree_format(form.getClass()) << std::endl;
			// disp(mat);
		}


		template<class Form>
		static void assemble_linear_and_print(const Form &form)
		{
			AssemblyContext<Backend> ctx;
			ctx.init_linear(form);

			ElementVector vec;

			FormEvaluator<Backend> eval;
			eval.eval(form, vec, ctx, true);

			// std::cout << tree_format(form.getClass()) << std::endl;
			// disp(vec);
		}

		static void run_interp_vec_test(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[run_interp_vec_test]" << std::endl;

			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
			auto V = Vx * Vy;

			auto u = trial(V);
			auto v = test(V);

			ElementVector c_uk = values(6, 1.);
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

			ElementVector c_uk = values(3, 1.);
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
		
		libMesh::MeshTools::Generation::build_square(*lm_mesh,
			10, 10,
			0, 1.,
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


		run_libmesh_test(init, [](
			LibMeshFormEvalTest &lm_test,
			const std::shared_ptr<libMesh::EquationSystems> &es) {
			
			//create system of equations
			auto &sys = es->add_system<libMesh::LinearImplicitSystem>("navier_stokes");

	
			auto Vx = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::SECOND, "vel_x");
			auto Vy = LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::SECOND, "vel_y");
			auto V  = Vx * Vy;

			auto Q =  LibMeshFunctionSpace(es, libMesh::LAGRANGE, libMesh::FIRST, "pressure");

			auto u = trial(V);
			auto v = test(V);

			auto ux = u[0];
			auto uy = u[1];

			auto p = trial(Q);
			auto q = test(Q);

			const double mu  = 1.;
			const double rho = 1.;

			DVectord sol;
			auto uk = interpolate(sol, u);
			auto g_uk = grad(uk);

			auto e = mu * (transpose(grad(u)) + grad(u));
			auto b_form_11 = integral(inner(e, grad(v)) + rho * inner(g_uk * u, v));
			auto b_form_12 = integral(inner(p, div(v)));
			auto b_form_21 = integral(inner(div(u), q));

			LMDenseVector r_v = zeros(2);
			auto l_form_1 =  integral(inner(coeff(0.), q));
			auto l_form_2 =  integral(inner(uk, v));

			auto lhs = b_form_11 + b_form_12 + b_form_21;
			auto rhs = l_form_1  + l_form_2;



			if(nl_solve(
				equations(
					lhs == rhs
				),
				constraints(
					boundary_conditions(ux == coeff(0.), {0, 1, 3}),
					boundary_conditions(ux == coeff(1.), {2}),
					boundary_conditions(uy == coeff(0.), {0, 1, 2, 3})
				),
				sol)) {
				//go back to libmesh
				convert(sol, *sys.solution);
				sys.solution->close();
				libMesh::ExodusII_IO(Q.mesh()).write_equation_systems ("navier_stokes.e", *es);
			} else {
				std::cerr << "[Error] solver failed to converge" << std::endl;
			}
		});

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
