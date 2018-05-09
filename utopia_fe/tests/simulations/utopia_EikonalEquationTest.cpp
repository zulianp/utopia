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
		const double diff_coeff = 1.;
		const double forcing_term = 5.;

		//discretization parameters
		const auto elem_type       = libMesh::QUAD4;
		const auto elem_order 	   = libMesh::FIRST;

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

		DVectord sol = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
		sol.set(1.);

		auto u_old = interpolate(sol, du);

		// Set (bi)linear forms
		auto l_form = (c1 * diff_coeff) * inner(inner(grad(u_old), grad(u_old)), v) * dX 
		               + c2 * inner(sqrt(inner(grad(u_old), grad(u_old))), v) * dX
		               - inner(coeff(forcing_term), v) * dX;

		auto b_form = (diff_coeff * c1) * inner(grad(du), grad(v)) * dX + c2 * inner(du * (1./sqrt(inner(u_old, u_old))), v) * dX;

		// assemble
		DSMatrixd hessian;
		DVectord gradient;

		utopia::assemble(b_form, hessian);
		utopia::assemble(l_form, gradient); 

		disp(hessian);
		disp(gradient);

	}
}
