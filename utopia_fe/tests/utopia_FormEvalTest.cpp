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

	class EqPreProcessor {
	public:
		EqPreProcessor()
		: n_bilinear_forms(0), n_linear_forms(0)
		{}

		template<class Any>
		inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }

		template<class Expr>
		inline int visit(Integral<Expr> &expr)
		{
			if(IsLinearForm<Expr>::value) {
				expr.set_integral_id(n_linear_forms++);
			} else if(IsBilinearForm<Expr>::value) {
				expr.set_integral_id(n_bilinear_forms++);
			}

			return TRAVERSE_CONTINUE;
		}

		template<class ExprTree>
		inline void apply(ExprTree &expr)
		{
			clear();
			traverse(expr, *this);
		}

		void clear()
		{
			n_linear_forms = 0;
			n_bilinear_forms = 0;
		}

	private:
		int n_bilinear_forms;
		int n_linear_forms;
	};

	template<int Backend>
	class ContextInitializer {
	public:
		template<class ExprTree>
		inline static void apply(const ExprTree &) {}
	};

	template<>
	class ContextInitializer<LIBMESH_TAG> {
	public:
		ContextInitializer(AssemblyContext<LIBMESH_TAG> &ctx)
		: ctx(ctx)
		{}

		template<class Any>
		inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }

		template<class Expr>
		inline int visit(Integral<Expr> &expr)
		{

		}

		template<class ExprTree>
		inline void apply(ExprTree &expr)
		{

			EqPreProcessor eq_pp;
			eq_pp.apply(expr);

			clear();
			traverse(expr, *this);
		}

		void clear()
		{
		}

	private:
		AssemblyContext<LIBMESH_TAG> &ctx;
	};



	class EqPrinter {
	public:
		template<class Eq>
		void operator()(const int index, const Eq &eq) const {
			std::cout << "equation: " << index << std::endl;
			std::cout << tree_format(eq.getClass()) << std::endl;
		}
	};





	//libmesh
	template<class FunctionSpaceT, class Left, class Right, class GlobalMatrix, class GlobalVector>
	void element_assemble_equation_v(const libMesh::MeshBase::const_element_iterator &it, const Equality<Left, Right> &equation, GlobalMatrix &mat, GlobalVector &vec)
	{
		typedef utopia::Traits<FunctionSpaceT> TraitsT;
		typedef typename TraitsT::Matrix ElementMatrix;
		typedef typename TraitsT::Vector ElementVector;

		static const int Backend = TraitsT::Backend;


		AssemblyContext<Backend> ctx;
		ctx.set_current_element((*it)->id());

		auto &&bilinear_form = equation.left();
		auto &&linear_form   = equation.right();

		ElementMatrix el_mat;
		ElementVector el_vec;

		ctx.init_bilinear(bilinear_form);

		FormEvaluator<Backend> eval;
		eval.eval(bilinear_form, el_mat, ctx, true);

		ctx.init_linear(linear_form);
		eval.eval(linear_form, el_vec, ctx, true);

		add_matrix(el_mat.implementation(), ctx.test_dof_indices, ctx.trial_dof_indices, mat);
		add_vector(el_vec.implementation(), ctx.test_dof_indices, vec);
	}

	//homemade
	template<class FunctionSpaceT, class Left, class Right, class Matrix, class Vector>
	void element_assemble_equation_v(const int element_index, const Equality<Left, Right> &equation, Matrix &mat, Vector &rhs)
	{
		std::cout << element_index << std::endl;
	}

	template<class FunctionSpaceT, class Left, class Right, class Matrix, class Vector>
	void assemble_equation_v(const Equality<Left, Right> &equation, Matrix &mat, Vector &vec)
	{
		auto &space = find_space<FunctionSpaceT>(equation);
		space.initialize();

		auto &m = space.mesh();
		auto &dof_map = space.dof_map();



		auto nnz_x_row = dof_map.n_local_dofs(); //FIXME
		mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
		vec = local_zeros(dof_map.n_local_dofs());

		Write<Matrix> w_m(mat);
		Write<Vector> w_v(vec);

		for(auto it = elements_begin(m); it != elements_end(m); ++it) {
			element_assemble_equation_v<FunctionSpaceT>(it, equation, mat, vec);
		}
	}

	template<class FunctionSpaceT, class Left, class Right>
	bool solve(const Equality<Left, Right> &equation, DVectord &sol)
	{
		DSMatrixd mat;
		DVectord vec;
		assemble_equation_v<FunctionSpaceT>(equation, mat, vec);
		sol = local_zeros(local_size(vec));

		//FIXME hardcoded boundary conditions
		{
			Write<DSMatrixd> w_(mat);
			mat.set(0, 0, 1.0);

			Size s = size(mat);
			for(uint i = 1; i < s.get(1); ++i) {
				mat.set(0, i, 0.0);
			}
		}

		disp(mat);
		disp(vec);

		Factorization<DSMatrixd, DVectord> solver;
		return solver.solve(mat, vec, sol);
	}


	class EqAssembler {
	public:
		EqAssembler(const libMesh::Elem * elem, DSMatrixd &mat, DVectord &vec) 
		: elem(elem), mat(mat), vec(vec)
		{
			
		}

		template<class Eq>
		void operator()(const int index, const Eq &eq_in) {
			typedef typename FindFunctionSpace<Eq>::Type FunctionSpaceT;
			typedef utopia::Traits<FunctionSpaceT> TraitsT;
			typedef typename TraitsT::Matrix ElementMatrix;
			typedef typename TraitsT::Vector ElementVector;

			Eq eq = eq_in;
			auto &space = find_space<FunctionSpaceT>(eq);

			auto &&bilinear_form = eq.left();
			auto &&linear_form   = eq.right();

			ElementMatrix el_mat;
			ElementVector el_vec;

			AssemblyContext<LIBMESH_TAG> ctx;

			ContextInitializer<LIBMESH_TAG> ci(ctx);
			ci.apply(eq);

			ctx.set_current_element(elem->id());
			ctx.init_bilinear(bilinear_form);

			FormEvaluator<LIBMESH_TAG> eval;
			eval.eval(bilinear_form, el_mat, ctx, true);

			ctx.init_linear(linear_form);
			eval.eval(linear_form, el_vec, ctx, true);

			std::vector<libMesh::dof_id_type> dof_indices;
			space.dof_map().dof_indices(elem, dof_indices);

			add_matrix(el_mat.implementation(), dof_indices, dof_indices, mat);
			add_vector(el_vec.implementation(), dof_indices, vec);
		};

		const libMesh::Elem * elem;
		DSMatrixd &mat;
		DVectord &vec;

		
	};


	template<class... Eqs, class... Constr>
	bool solve(const Equations<Eqs...> &eqs, const FEConstraints<Constr...> &constr, DVectord &sol)
	{

		//FIXME this stuff only works for the libmesh backend
		typedef typename GetFirst<Eqs...>::Type Eq1Type;
		typedef typename FindFunctionSpace<Eq1Type>::Type FunctionSpaceT;
		auto &space = find_space<FunctionSpaceT>(eqs.template get<0>());


		
		space.initialize();
		auto &m = space.mesh();
		auto &dof_map = space.dof_map();
		auto nnz_x_row = dof_map.n_local_dofs(); //FIXME


		DSMatrixd mat;
		DVectord vec;

		mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
		vec = local_zeros(dof_map.n_local_dofs());

		{
			Write<DSMatrixd> w_m(mat);
			Write<DVectord>  w_v(vec);

			for(auto it = elements_begin(m); it != elements_end(m); ++it) {
				EqAssembler eq_assembler(*it, mat, vec);
				eqs.each(eq_assembler);
			}
		}	


		disp(mat);
		disp(vec);

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
			EqPreProcessor preproc;

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

			assemble_bilinear_and_print(b_form_11 + b_form_12 + b_form_21);
			assemble_linear_and_print(l_form_1 + l_form_2);

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

			ElementMatrix A = identity(2, 2);
			{
				Write<ElementMatrix> w(A);
				A.set(1, 1, 2.);
			}

			auto mixed = integral( inner(grad(u), 0.1 * A * v) );
			assemble_bilinear_and_print(mixed);
		}

		static void leastsquares_helmoholtz(const std::shared_ptr<SpaceInput> &space_input)
		{
			std::cout << "[leastsquares_helmoholtz]" << std::endl;

			//scalar
			auto U = FunctionSpaceT(space_input);

			//vector
			auto Vx = FunctionSpaceT(space_input);
			auto Vy = FunctionSpaceT(space_input);
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
			// assemble_equation_v<FunctionSpaceT>(eq, mat, vec);
		}
	};

	typedef FormEvalTest<libMesh::EquationSystems, utopia::LibMeshFunctionSpace> LibMeshFormEvalTest;

	template<typename F>
	static void run_libmesh_test(libMesh::LibMeshInit &init, F fun)
	{	
		auto lm_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		
		libMesh::MeshTools::Generation::build_square(*lm_mesh,
			1, 1,
			0, 1.,
			0, 1.,
			libMesh::TRI3);

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

		run_libmesh_test(init,[](
			LibMeshFormEvalTest &lm_test,
			const std::shared_ptr<libMesh::EquationSystems> &es) {
			es->add_system<libMesh::LinearImplicitSystem>("run_navier_stokes_test");
			lm_test.run_navier_stokes_test(es);
		});

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
			
		// 	//create system of equations
		// 	es->add_system<libMesh::LinearImplicitSystem>("run_local_2_global_test");

		// 	auto V = LibMeshFunctionSpace(es);

		// 	auto u = trial(V);
		// 	auto v = test(V);

		// 	// DVectord sol;
		// 	// bool success = solve<LibMeshFunctionSpace>(
		// 	// 	inner(grad(u), grad(v)) * dX == inner(coeff(1.), v) * dX,
		// 	// 	sol);

		// 	/*
		// 		bool success = solve(
		// 			equations(
		// 				inner(grad(u), grad(v)) * dX == inner(coeff(1.), v) * dX,
		// 				...
		// 			),
		// 			constraints(
		// 				//dirichlet
		// 				boundary_conditions(u == 0., {1}),
		// 				boundary_conditions(u == 1., {2}),
		// 				//neumann (?)
		// 				boundary_conditions(inner(grad(u), normal_) == x_ * y_, {3}),
		// 				boundary_conditions(inner(grad(u), normal_) == inner(func, normal_), {3}),
		// 				//pseudo l2-projection
		// 				inner(u, biorth(v)) * dV({1}) == inner(w, biorth(v)) * dV({2})
		// 				//or
		// 				inner(u - w, biorth(v)) * dX == 0,
		// 				//contact
		// 				contact_conditions({{3, 4}}, search_radius)

		// 			),
		// 			sol
		// 		);
		// 	*/
			

		// 	// std::cout << "solved: " << (success ? "true" : "false") << std::endl;
		// 	// disp(sol);
		// });

		// run_libmesh_test(init, [](
		// 	LibMeshFormEvalTest &lm_test,
		// 	const std::shared_ptr<libMesh::EquationSystems> &es) {
			
		// 	//create system of equations
		// 	es->add_system<libMesh::LinearImplicitSystem>("test_equations");

		// 	//space for u
		// 	auto V = LibMeshFunctionSpace(es);

		// 	auto u = trial(V);
		// 	auto v = test(V);

		// 	//space for gradient of u
		// 	auto W1 = LibMeshFunctionSpace(es);
		// 	auto W2 = LibMeshFunctionSpace(es);

		// 	auto W = W1 * W2;

		// 	auto q = trial(W);
		// 	auto w = test(W);

		// 	LMDenseVector rhs_w = values(2, 1.);

		// 	DVectord sol;
		// 	const bool success = solve(
		// 		equations(
		// 			inner(grad(u), grad(v)) * dX == inner(coeff(1.), v) * dX,
		// 			// inner(grad(u), w) * dX       == inner(coeff(rhs_w), v) * dX,
		// 			inner(q, w) * dX 			 == inner(coeff(rhs_w), w) * dX
		// 		),
		// 		constraints(
		// 			boundary_conditions(u == coeff(-0.2), {1}),
		// 			boundary_conditions(u == coeff(0.2),  {2})
		// 		),
		// 		sol);
		// });


		/////////////////////////////////////////////////////////////////////
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
