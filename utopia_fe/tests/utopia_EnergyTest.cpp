#include "utopia_EnergyTest.hpp"

#include "utopia_libmesh.hpp"


namespace utopia {
	void run_energy_test(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());


		const int n = 1;
		libMesh::MeshTools::Generation::build_square(*mesh,
			n, n,
			0, 1,
			0, 1.,
			// libMesh::QUAD4
			libMesh::TRI3
			);

		auto dim = mesh->mesh_dimension();

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("eval-test");

		const double mu = 0.1;
		const double lambda = 0.1;

		auto elem_order = libMesh::FIRST;

		////////////////////////////////////////////

		auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
		auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
		auto V = Vx * Vy;

		auto u = trial(V);
		auto v = test(V);

		Vx.initialize();
		DVectord sol = ghosted(Vx.dof_map().n_local_dofs(), Vx.dof_map().n_dofs(), Vx.dof_map().get_send_list());

		// auto uh = interpolate(V); 

	}
}