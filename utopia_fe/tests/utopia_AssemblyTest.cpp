#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"

namespace utopia {
	void run_assembly_test(libMesh::LibMeshInit &init)
	{
		std::cout << "[run_assembly_test]" << std::endl;
		typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());

		const unsigned int n = 2;
		libMesh::MeshTools::Generation::build_cube(*mesh,
			n, n, n,
			0, 1,
			0, 1.,
			0, 1.,
			libMesh::TET4);

		auto es = std::make_shared<libMesh::EquationSystems>(*mesh);
		es->add_system<libMesh::LinearImplicitSystem>("lapl");

		auto V = FunctionSpaceT(es);

		auto u = trial(V);
		auto v = test(V);


		//set-up boundary conditions
		//...
		//... add boundary conditions to system

		V.initialize();

		const double alpha = 1.;
		USparseMatrix laplacian, mass_matrix;

		assemble(inner(alpha * grad(u), grad(v)) * dX, laplacian);
		assemble(inner(u, v) * dX, mass_matrix);

		const double norm_laplacian   = norm2(laplacian);
		const double norm_mass_matrix = norm2(mass_matrix);

		assert(approxeq(9.60035,   norm_laplacian, 1e-5));
		assert(approxeq(0.0660453, norm_mass_matrix, 1e-5));
	}
}
