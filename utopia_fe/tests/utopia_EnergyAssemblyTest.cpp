#include "utopia_EnergyAssemblyTest.hpp"
#include "utopia_libmesh.hpp"

#include <algorithm>

namespace utopia {

	static void laplacian_test(LibMeshFunctionSpace &V)
	{
		UVector x = local_values(V.dof_map().n_local_dofs(), 10.);
		auto rhs = coeff(1.);

		auto u  = trial(V);
		auto v  = test(V);

		auto uk = interpolate(x, u);

		auto f1 = 0.5 * inner(grad(uk), grad(uk)) * dX;
		auto f2 = inner(rhs, uk) * dX;

		double energy1 = -1.;
		utopia::assemble(f1, energy1);

		double energy2 = -1.;
		utopia::assemble(f2, energy2);

		double energy = energy1 - energy2;

		utopia_test_assert(std::abs(energy1) < 1e-8);

		std::cout << std::abs(energy) << std::endl;

		USparseMatrix H;
		UVector g;

		utopia::assemble(inner(grad(u), grad(v)) * dX == inner(rhs, v) * dX, H, g);
		apply_boundary_conditions(V.dof_map(), H, g);

		Factorization<USparseMatrix, UVector> solver;
		solver.solve(H, g, x);

		double energy_after = -1.;
		utopia::assemble(
			// inner_volume_only_integral(
			integral(
				0.5 * inner(grad(uk), grad(uk)) - inner(rhs, uk)
			), 
			energy_after
		);

		utopia_test_assert(std::abs(energy) > std::abs(energy_after));
		std::cout << std::abs(energy_after) << std::endl;

		write("test.e", V, x);
	}

	void run_energy_test(libMesh::LibMeshInit &init)
	{
		libMesh::DistributedMesh mesh(init.comm());

		int n = 10;
		libMesh::MeshTools::Generation::build_square(
			mesh,
			n,  n,
			0., 1.,
			0., 1.,
			libMesh::QUAD4
		);

		LibMeshFunctionSpace V(mesh); 

		auto u = trial(V);
		init_constraints( constraints(
			boundary_conditions(u == coeff(0), {1, 3}),
			// boundary_conditions(u == coeff(1), {0, 2})
			boundary_conditions(u == coeff(0), {0, 2})
		));
		

		V.initialize();


		//////////////////////////////////////////
		//////////////////////////////////////////

		laplacian_test(V);

		//////////////////////////////////////////
	}
}
