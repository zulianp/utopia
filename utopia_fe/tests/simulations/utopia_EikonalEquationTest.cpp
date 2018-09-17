#include "utopia_LeastSquaresHelmholtz.hpp"

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

namespace utopia {

	void run_eikonal_equation_test(libMesh::LibMeshInit &init)
	{
		typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

		//model parameters
		const unsigned int n = 2;
		const double c1 = 1.;
		const double c2 = 1.;
		const double tau = 1.;
		const double diffusivity = 1.;
		const double forcing_term = 5.;

#ifdef WITH_TINY_EXPR
		auto f = symbolic("5 * sqrt(x^2 + y^2 + z^2)");
#else
		auto f = coeff(forcing_term);
#endif  //WITH_TINY_EXPR

		//discretization parameters
		const auto elem_type = libMesh::QUAD8;
		const auto elem_order = libMesh::SECOND;

		//mesh
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
		libMesh::MeshTools::Generation::build_square(
			*mesh,
			n, n,
			0, 1,
			0, 1.,
			elem_type
			);

		//equations system
		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("eikonal_equation");

		//scalar function space
		auto V = FunctionSpaceT(equation_systems, libMesh::LAGRANGE, elem_order, "u");


		auto du = trial(V);
		auto v  = test(V);

		V.initialize();
		auto &dof_map = V.dof_map();
		dof_map.prepare_send_list();

		UVector sol = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
		sol.set(0.);

		UVector diff_coeff = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
		diff_coeff.set(diffusivity);

		auto u_old = interpolate(sol, du);
		auto d     = interpolate(diff_coeff, du);

		// Set (bi)linear forms
		auto l_form = c1 * inner(d * grad(u_old), grad(v)) * dX
		+ c2 * inner(sqrt(inner(grad(u_old), grad(u_old))), v) * dX
		- inner(f, v) * dX;

		auto b_form = c1 * inner(d * grad(du), grad(v)) * dX
		+ c2 * inner(inner(grad(du), grad(u_old))/(coeff(1e-10) + sqrt(inner(grad(u_old), grad(u_old)))), v) * dX;

		// assemble
		USparseMatrix hessian;
		UVector gradient;

		utopia::assemble(b_form, hessian);
		utopia::assemble(l_form, gradient);

		disp(hessian);
		disp(gradient);

	}
}
