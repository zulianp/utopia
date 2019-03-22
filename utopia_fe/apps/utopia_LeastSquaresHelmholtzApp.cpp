#include "utopia_libmesh.hpp"
#include "utopia_LeastSquaresHelmholtzApp.hpp"

#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

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

namespace utopia {

	void LeastSquaresHelmholtzApp::run(Input &in)
	{
		typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

		//model parameters
		const unsigned int n = 50;
		const double c = 10.0;
		const double beta = 0.99;
		const double f = 5.;

		//discretization parameters
		const auto elem_type       = libMesh::QUAD8;
		const auto elem_order 	   = libMesh::SECOND;
		const auto elem_order_grad = libMesh::FIRST;

		//mesh
		auto mesh = std::make_shared<libMesh::DistributedMesh>(comm());
		libMesh::MeshTools::Generation::build_square(
			*mesh,
			n, n,
			0, 1,
			0, 1.,
			elem_type
		);

		//equations system
		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("leastsquares_helmoholtz");

		//scalar function space
		auto V = FunctionSpaceT(equation_systems, libMesh::LAGRANGE, elem_order, "u");

		//vector function space
		auto Qx = FunctionSpaceT(equation_systems, libMesh::LAGRANGE, elem_order_grad, "grad_x");
		auto Qy = FunctionSpaceT(equation_systems, libMesh::LAGRANGE, elem_order_grad, "grad_y");
		auto Q  = Qx * Qy;

		auto u = trial(V);
		auto v = test(V);

		auto s = trial(Q);
		auto q = test(Q);

		auto sx = s[0];
		auto sy = s[1];

		//bilinear forms
		auto b_11 = integral((c*c) * inner(u, v) + inner(grad(u), grad(v)));
		auto b_12 = integral(c * inner(div(s), v) + inner(s, grad(v)));
		auto b_21 = integral(c * inner(u, div(q)) + inner(grad(u), q));
		auto b_22 = integral(inner(s, q) + inner(div(s), div(q)) + beta * inner(curl(s), curl(q)));

		//linear forms
		auto l_1 = integral(c * inner(coeff(f), v));
		auto l_2 = integral(inner(coeff(f), div(q)));

		auto b_form = b_11 + b_12 + b_21 + b_22;
		auto l_form = l_1 + l_2;

		UVector sol;
		if(!solve(
			equations(
				b_form == l_form
			),
			constraints(
				boundary_conditions(u  == coeff(0.0), {0, 1, 2}),
				boundary_conditions(sx == coeff(0.0), {3}),
				boundary_conditions(sy == coeff(0.0), {3})
			),
			sol)) {

			std::cerr << "[Error] unable to solve" << std::endl;
			return;
		}


		libMesh::ExodusII_IO io(*mesh);
		convert(sol, *sys.solution);
		io.write_equation_systems("leastsquares_helmoholtz.e", *equation_systems);
	}
}
