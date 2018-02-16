#include "utopia_SteadyContact.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_Newmark.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"

namespace utopia {
	template class ContactSolver<DSMatrixd, DVectord>;

	// typedef utopia::ContactSolver<DSMatrixd, DVectord> ContactSolverT;
	// typedef utopia::Newmark<DSMatrixd, DVectord> ContactSolverT;
	typedef utopia::ContactStabilizedNewmark<DSMatrixd, DVectord> ContactSolverT;

	void run_steady_contact(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
		mesh->read("../data/wear_2_far.e");

		const auto dim = mesh->mesh_dimension();

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);	
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("dynamic-contact");

		const double dt = 0.1;
		LameeParameters lamee_params(20., 20.);
		// lamee_params.set_mu(2, 50.);
		// lamee_params.set_lambda(2, 100.);


		auto elem_order = libMesh::FIRST;

		////////////////////////////////////////////

		auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
		auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
		auto V = Vx * Vy;

		// if(dim == 3) {
		// 	V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
		// }


		auto u = trial(V);
		auto ux = u[0];
		auto uy = u[1];


		auto constr = constraints(
			boundary_conditions(ux == coeff(0.),    {4}),
			boundary_conditions(uy == coeff(0.),    {4})
		);

		init_constraints(constr);
		Vx.initialize();


		auto ef = std::make_shared<ConstantExternalForce>();

		auto vx = test(Vx);
		auto vy = test(Vy);

		ef->init(integral(inner(coeff(0.), vx) + inner(coeff(-.2), vy), 1));
		// ef->init(integral(inner(coeff(0.), vx) + inner(coeff(-.2), vy)));

		auto material = std::make_shared<NeoHookean<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);
		// auto material = std::make_shared<IncompressibleNeoHookean<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);
		// auto material = std::make_shared<SaintVenantKirchoff<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);
		// auto material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);

		ContactParams contact_params;
		// contact_params.contact_pair_tags = {{2, 1}};
		contact_params.contact_pair_tags = {{1, 2}};
		contact_params.search_radius = 0.3;

		ContactSolverT sc(make_ref(V), material, dt, contact_params); 
		sc.set_external_force_fun(ef);		
		sc.initial_condition();
		sc.solve_dynamic(400);
	}
}
