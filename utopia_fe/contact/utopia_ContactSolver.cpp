#include "utopia_ContactSolver.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_Newmark.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"


#include "libmesh/mesh_refinement.h"

namespace utopia {
	template class ContactSolver<DSMatrixd, DVectord>;
	typedef utopia::ProductFunctionSpace<LibMeshFunctionSpace> VectorFunctionSpace;

	// typedef utopia::ContactSolver<DSMatrixd, DVectord> ContactSolverT;
	// typedef utopia::Newmark<DSMatrixd, DVectord> ContactSolverT;
	typedef utopia::ContactStabilizedNewmark<DSMatrixd, DVectord> ContactSolverT;

	void run_steady_contact(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
		// mesh->read("../data/wear_2_far.e");
		// mesh->read("../data/channel_2d.e");
		mesh->read("../data/leaves_3d.e");

		// {
		// 	libMesh::MeshRefinement mesh_refinement(*mesh);
		// 	mesh_refinement.make_flags_parallel_consistent();
		// 	mesh_refinement.uniformly_refine(2);
		// }

		const auto dim = mesh->mesh_dimension();

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);	
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("dynamic-contact");

		double dt = 0.05;
		if(dim == 3) {
			dt = 0.001;
		}
		
		// LameeParameters lamee_params(20., 20.);
		// lamee_params.set_mu(2, 10.);
		// lamee_params.set_lambda(2, 10.);

		LameeParameters lamee_params(1., 1.);

		auto elem_order = libMesh::FIRST;

		////////////////////////////////////////////

		auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
		auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
		auto V = Vx * Vy;

		// auto W = VectorFunctionSpace(dim, equation_systems, libMesh::LAGRANGE, elem_order);

		if(dim == 3) {
			V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
		}

		auto u = trial(V);
		auto ux = u[0];
		auto uy = u[1];

		auto constr = constraints(
			boundary_conditions(ux == coeff(0.), {4}),
			boundary_conditions(uy == coeff(0.), {4})
		);

		if(dim == 3) {
			auto uz = u[2];
			auto constr3 = constr + boundary_conditions(uz == coeff(0.), {4});
			init_constraints(constr3);
		} else {
			init_constraints(constr);
		}

		Vx.initialize();

		auto ef = std::make_shared<ConstantExternalForce>();

		auto vx = test(Vx);
		auto vy = test(Vy);

		// ef->init(integral(inner(coeff(0.), vx) + inner(coeff(-.2), vy), 1));
		
		if(dim == 3) {
			ef->init(integral(inner(coeff(7.), vx)));
		} else {
			ef->init(integral(inner(coeff(-.2), vy)));	
		}

		// auto material = std::make_shared<NeoHookean<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);
		// auto material = std::make_shared<IncompressibleNeoHookean<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);
		// auto material = std::make_shared<SaintVenantKirchoff<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);
		auto material = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, lamee_params);

		ContactParams contact_params;
		// contact_params.contact_pair_tags = {{2, 1}};
		contact_params.contact_pair_tags = {{1, 2}, {1, 3}, {2, 3}};

		if(dim == 3) {
			contact_params.search_radius = 0.0001;
		} else {
			contact_params.search_radius = 0.1;
		}

		
		ContactSolverT sc(make_ref(V), material, dt, contact_params); 
		sc.set_tol(5e-6);
		
		//begin: multigrid
		auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
		auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord> >();
		auto mg = std::make_shared<SemiGeometricMultigrid>(smoother, linear_solver);
		mg->verbose(true);
		mg->init(Vx, 4);
		
		mg->algebraic().atol(1e-15);
		mg->algebraic().rtol(1e-15);
		mg->algebraic().stol(1e-13);

		sc.set_linear_solver(mg);
		//end: multigrid



		sc.set_external_force_fun(ef);		
		sc.initial_condition();
		sc.solve_dynamic(400);
	}
}
