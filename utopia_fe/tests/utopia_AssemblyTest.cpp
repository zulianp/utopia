#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_partitioner.h"

namespace utopia {
	void run_assembly_test(libMesh::LibMeshInit &init)
	{
		std::cout << "[run_assembly_test]" << std::endl;
		typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());	
		// const std::string data_path = Utopia::Instance().get("data_path");
		// const std::string mesh_path = data_path + "/fine_mesh.e";
		// mesh->read(mesh_path);

		const unsigned int n = 2;
		libMesh::MeshTools::Generation::build_cube(*mesh,
			n, n, n,
			0, 1,
			0, 1.,
			0, 1.,
			libMesh::TET4);


		auto es = std::make_shared<libMesh::EquationSystems>(*mesh);
		auto &sys =es->add_system<libMesh::LinearImplicitSystem>("lapl");
		auto V = FunctionSpaceT(es);
		V.initialize();

		auto u = trial(V);
		auto v = test(V);

		const double alpha = 1.;
		auto laplacian = inner(alpha * grad(u), grad(v)) * dX;

		DSMatrixd mat, lm_mat;
		assemble(laplacian, mat);
		assemble(laplacian, *sys.matrix);
		sys.matrix->close();
		convert(*sys.matrix, lm_mat);
		
		double sum_lm_mat = norm2(lm_mat);
		double sum_mat = norm2(mat);
		std::cout << "sum_mat: " << sum_lm_mat << " == " << sum_mat << std::endl;
	}
}