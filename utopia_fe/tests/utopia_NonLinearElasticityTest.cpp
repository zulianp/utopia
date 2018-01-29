#include "utopia_NonLinearElasticityTest.hpp"

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
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"

#include "libmesh/exodusII_io.h"
#include <algorithm>

namespace utopia {
	template<class Expr>
	Binary<Expr, Number<double>, Minus> operator-(const Expression<Expr> &left, const double &right)
	{
		return Binary<Expr, Number<double>, Minus>(left, right); 
	}

	template<class Expr>
	double operator-(const Number<double> &left, const Number<double> &right)
	{
		return static_cast<double>(left) - static_cast<double>(right);
	}

	void run_non_linear_elasticity_test(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());

		libMesh::MeshTools::Generation::build_square(*mesh,
			10, 10,
			0, 1,
			0, 1.,
			libMesh::QUAD8);

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);	
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("neo-hookean");

		const double mu = 1.;
		const double lambda = 1.;

		////////////////////////////////////////////

		auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::SECOND, "disp_x");
		auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::SECOND, "disp_y");
		auto V = Vx * Vy;

		auto u = trial(V);
		auto v = test(V);

		auto ux = u[0];
		auto uy = u[1];

		////////////////////////////////////////////

		DVectord sol;
		auto uk = interpolate(sol, u);

		auto F 		 = identity() + grad(uk);
		auto F_inv   = inv(F);
		auto F_inv_t = transpose(F_inv);
		auto g_uk    = grad(uk);

		auto J = det(F);
		auto C = F + transpose(F);
		auto t = trace(C);
		
		auto P = mu * (F - F_inv_t) + lambda * logn(J) * F_inv_t;
		
		//compressible neo-hookean
		auto l_form = inner(P, grad(v)) * dX;
		auto b_form = (
			mu * inner(grad(u), grad(v))
			- inner((lambda * logn(J) - mu) * transpose(F_inv * grad(u)), F_inv * grad(v))
			+ inner(lambda * F_inv_t, grad(u)) * inner(F_inv_t, grad(v))
			) * dX;

		////////////////////////////////////////////

		if(nl_solve(
			equations(
				b_form == l_form
				),
			constraints(
				boundary_conditions(ux == coeff(0.), {1, 3}),
				boundary_conditions(uy == coeff(-0.1), {1}),
				boundary_conditions(uy == coeff(0.1),  {3})
				),
			sol)) {

			convert(sol, *sys.solution);
			sys.solution->close();
			libMesh::ExodusII_IO(*mesh).write_equation_systems("neo-hookean.e", *equation_systems);
	} else {
		std::cerr << "[Error] solver failed to converge" << std::endl;
	}

}

}