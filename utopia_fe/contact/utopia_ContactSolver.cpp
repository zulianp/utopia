#include "utopia_ContactSolver.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_Newmark.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"


#include "libmesh/mesh_refinement.h"

namespace utopia {
	template class ContactSolver<USMatrix, UVector>;
	typedef utopia::ProductFunctionSpace<LibMeshFunctionSpace> VectorFunctionSpace;

	// typedef utopia::ContactSolver<USMatrix, UVector> ContactSolverT;
	// typedef utopia::Newmark<USMatrix, UVector> ContactSolverT;
	typedef utopia::ContactStabilizedNewmark<USMatrix, UVector> ContactSolverT;

	void run_steady_contact(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
		// mesh->read("../data/wear_2_b.e");
		mesh->read("../data/wear_2_far.e");
		// mesh->read(utopia::Utopia::instance().get("data_path") + "/input_file.e");
		// mesh->read("../data/channel_2d.e");
		// mesh->read("../data/leaves_3d_b.e");


		const auto dim = mesh->mesh_dimension();

		// if(dim == 2)
		// {
			// libMesh::MeshRefinement mesh_refinement(*mesh);
			// mesh_refinement.make_flags_parallel_consistent();
			// mesh_refinement.uniformly_refine(2);
		// }

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("dynamic-contact");

		double dt = 0.05;
		if(dim == 3) {
			dt = 0.0001;
		}

		double mu = 20., lambda = 20.;
		if(dim == 3) {
			mu = 300.;
			lambda = 600.;
		}

		LameeParameters lamee_params(mu, lambda);

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

		std::cout << "n_dofs: " << Vx.dof_map().n_dofs() << std::endl;

		auto ef = std::make_shared<ConstantExternalForce>();

		auto vx = test(Vx);
		auto vy = test(Vy);

		// ef->init(integral(inner(coeff(0.), vx) + inner(coeff(-.2), vy), 1));

		if(dim == 3) {
			ef->init(integral(inner(coeff(7000.), vx)));
		} else {
			ef->init(integral(inner(coeff(-.2), vy)));
		}

		// auto material = std::make_shared<NeoHookean<decltype(V), USMatrix, UVector>>(V, lamee_params);
		// auto material = std::make_shared<IncompressibleNeoHookean<decltype(V), USMatrix, UVector>>(V, lamee_params);
		// auto material = std::make_shared<SaintVenantKirchoff<decltype(V), USMatrix, UVector>>(V, lamee_params);
		auto material = std::make_shared<LinearElasticity<decltype(V), USMatrix, UVector>>(V, lamee_params);

		ContactParams contact_params;
		// contact_params.contact_pair_tags = {{2, 1}};
		contact_params.contact_pair_tags = {{1, 2}, {1, 3}, {2, 3}};

		if(dim == 3) {
			contact_params.search_radius = 0.0001;
		} else {
			contact_params.search_radius = 0.04;
		}


		// auto stabilized_material = std::make_shared<StabilizedMaterial<decltype(V), USMatrix, UVector> >(V, 1e-2, material);
		// ContactSolverT sc(make_ref(V), stabilized_material, dt, contact_params);

		ContactSolverT sc(make_ref(V), material, dt, contact_params);
		sc.set_tol(5e-6);

		// auto ls = std::make_shared<Factorization<USMatrix, UVector>>();
		// auto ls = std::make_shared<GMRES<USMatrix, UVector>>();
		// ls->atol(1e-15);
		// ls->rtol(1e-15);
		// ls->stol(1e-15);
		// ls->max_it(1000);
		// // ls->verbose(true);
		// sc.set_linear_solver(ls);
		// sc.set_bypass_contact(true);
		sc.set_max_outer_loops(30);

		// begin: multigrid

		// auto linear_solver = std::make_shared<BiCGStab<USMatrix, UVector>>();
		// auto smoother = std::make_shared<ConjugateGradient<USMatrix, UVector, HOMEMADE>>();
		// auto smoother = std::make_shared<BiCGStab<USMatrix, UVector>>();
		// prec->max_it(1);
		// smoother->set_preconditioner(prec);

		// auto smoother = std::make_shared<SOR<USMatrix, UVector> >();
		// // auto smoother = std::make_shared<GMRES<USMatrix, UVector> >();

		// auto linear_solver = std::make_shared<Factorization<USMatrix, UVector>>();
		// // auto linear_solver = std::make_shared<BiCGStab<USMatrix, UVector>>();
		// // auto smoother = std::make_shared<ProjectedGaussSeidel<USMatrix, UVector, HOMEMADE> >();
		// auto mg = std::make_shared<SemiGeometricMultigrid>(smoother, linear_solver);
		// mg->verbose(true);
		// mg->set_use_interpolation(true);
		// mg->init(Vx, 4);


		// mg->algebraic().atol(1e-18);
		// mg->algebraic().rtol(1e-8);
		// mg->algebraic().stol(1e-16);
		// // mg->algebraic().set_use_line_search(true);

		// sc.set_linear_solver(mg);

		// mg->set_separate_subdomains(true);
		// end: multigrid

		sc.set_external_force_fun(ef);
		sc.initial_condition(2.);
		sc.solve_dynamic(400);
	}
}
