
#include <iostream>
#include <cmath>
#include <iomanip>
#include "utopia_FEDSLBaseExamples.hpp"
#include "utopia.hpp"

//fe extension
#include "utopia_TimeDiffExamples.hpp"
#include "utopia_fe_core.hpp"
#include "MortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>

using namespace libMesh;
using namespace std;

namespace utopia {

#ifndef WITHOUT_TIME_DIFF_EXAMPLES	

	template<class Expr_>
	class ExplicitIntegratorTransform {
	public:
		typedef Expr_ Expr;

		const Expr &apply(const Expr &expr)
		{
			return expr;
		}
	};

	template<class Left, class Right>
	class ExplicitIntegratorTransform< Equality<Left, Right> > {
	public:
		typedef Equality<Left, Right>  Expr;
		typedef typename Left::Expr::Expr::Right TestFun;
		typedef typename Left::Expr::Expr::Left::Expr Solution;
		typedef utopia::Reduce<Binary<TestFun, TestFun, EMultiplies>, Plus> L2Inner;
		typedef utopia::Multiply<Inverse<Integral<L2Inner>>, Right> UpdateFunction;

		inline TestFun &get_test_fun(const Expr &expr)
		{
			return const_cast<TestFun &>(expr.left().expr().expr().right());
		}

		inline Solution &get_solution(const Expr &expr)
		{
			return const_cast<Solution &>(expr.left().expr().expr().left().expr());
		}

		inline UpdateFunction get_update_function(const Expr &expr)
		{
			auto &test_fun = get_test_fun(expr);
			return inv( integral(dot(test_fun, test_fun)) ) * expr.right();
		}
	};

	template<class Left, class Right, class Context>
	void integrate(const Equality<Left, Right> &eq, const double t0, const double delta_t, const double t_end, Context &context)
	{
		ExplicitIntegratorTransform< Equality<Left, Right> > trafo;
		auto &v = trafo.get_test_fun(eq);
		auto uf = trafo.get_update_function(eq);

		ExodusII_IO io(v.mesh());
		
		auto mass 	= uf.left().expr();
		auto l_form = uf.right();

		auto s = size(v);
		DVectord g  = zeros(s);
		DSMatrixd m = sparse(s.get(0), s.get(0), 0.2 * s.get(0));
		DSMatrixd boundary_matrix = identity(s.get(0), s.get(0));

		{	
			Write<DSMatrixd> w_m(m);
			Write<DVectord>  w_g(g);
			assemble(v, v, mass, m, false);
			assemble(v, l_form,  g, false);
		}

		DVectord Minvg = zeros(s);
		solve(m, g, Minvg);

		auto &sol = trafo.get_solution(eq);
		DVectord u_sol = zeros(s);
		convert(sol.vector(), u_sol);
		apply_boundary_conditions(v, boundary_matrix, u_sol);
		u_sol += delta_t * Minvg;
		apply_boundary_conditions(v, boundary_matrix, u_sol);

		double t_range = t_end - t0;
		uint T = (t_range)/delta_t;

		int write_iteration = 1;
		for(uint t = 1; t < T; ++t) {

			convert(u_sol, sol.vector());
			if(t % 20 == 0) {
				std::stringstream ss;
			 	io.write_timestep("heat_equation.e", context.equation_systems, write_iteration++, t);
			}
			{	
				Write<DVectord> w_g(g);
				assemble(v, l_form, g);
			}

			solve(m, g, Minvg);
			u_sol += delta_t * Minvg;
		}
	}

	void run_time_diff_examples(libMesh::LibMeshInit &init)
	{
		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
		MeshTools::Generation::build_square (*mesh,
			10, 10,
			-1., 1.,
			-1., 1.,
			QUAD4);

		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		auto Vh = fe_space(LAGRANGE, FIRST, context);
		auto v  = fe_function(Vh);

		strong_enforce( boundary_conditions(v == coeff(.0), {1, 2}) );


		const int dim = mesh->mesh_dimension();
		std::function<Real(const Point &p) > rhs_fun = [dim](const Point &p) -> Real {
			Real center = 0;

			Real ret = 0;
			for(int i = 0; i < dim; ++i) {
				ret += (p(i) - center) * (p(i) - center);
			}

			return std::sqrt(ret);
		};

		context.equation_systems.init();

		v.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));

		auto v_k = interpolate(coeff(.5), v, make_ref(*context.system.solution));
		auto eq  = integral(dot(dt(v_k), v)) == integral(dot(coeff(rhs_fun), v) -  0.1 * dot(grad(v_k), grad(v)));

		auto ass = make_assembly([&]() -> void {
			assemble(v, v, integral(dot(v, v)), integral(dot(coeff(0.0), v)), *context.system.matrix, *context.system.rhs);
		});

		context.system.attach_assemble_object(ass);
		context.equation_systems.print_info();
		context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 0;
		context.equation_systems.solve();
		integrate(eq, 0, 0.001, 2, context); 		
	}

#else
	void run_time_diff_examples(libMesh::LibMeshInit &init) {}
#endif	

}