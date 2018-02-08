#include "utopia_SteadyContact.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

namespace utopia {
	template class SteadyContact<DSMatrixd, DVectord>;
	typedef utopia::SteadyContact<DSMatrixd, DVectord> SteadyContactT;


	void run_steady_contact(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
		mesh->read("../data/wear_2_far.e");

		const auto dim = mesh->mesh_dimension();

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);	
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("steady-contact");

		const double mu = 1.;
		const double lambda = 1.;

		auto elem_order = libMesh::FIRST;

		////////////////////////////////////////////

		auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "u_x");
		auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "u_y");
		auto V = Vx * Vy;

		// if(dim == 3) {
		// 	V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
		// }


		auto u = trial(V);
		auto ux = u[0];
		auto uy = u[1];


		auto constr = constraints(
			boundary_conditions(uy == coeff(0.32),  {4}),
			boundary_conditions(uy == coeff(-0.28), {3}),
			boundary_conditions(ux == coeff(0.),    {3, 4})
		);

		init_constraints(constr);
		Vx.initialize();


		auto material = std::make_shared<NeoHookean<decltype(V), DSMatrixd, DVectord>>(V, mu, lambda);
		ContactParams contact_params;
		contact_params.contact_pair_tags = {{1, 2}};
		contact_params.search_radius = 1.7;

		SteadyContactT sc(make_ref(V), material, contact_params);
		sc.solve();

	}
}
