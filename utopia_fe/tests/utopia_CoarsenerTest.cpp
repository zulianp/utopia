#include "utopia_CoarsenerTest.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_BoundingBoxCoarsener.hpp"
#include "utopia_Socket.hpp"

#include "libmesh/exodusII_io.h"

namespace utopia {
	void run_coarsener_test(libMesh::LibMeshInit &init)
	{
		std::cout << "[run_coarsener_test]" << std::endl;
		typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());	
		// mesh->read("../data/leaves_3d_b.e");
		mesh->read("../data/wear_2_far.e");

		if(mpi_world_size() == 1) {
			plot_mesh(*mesh, "mesh/fine");
		}

		BoundingBoxCoarsener bb_coarsener;
		bb_coarsener.init(4, *mesh);
		bb_coarsener.describe();

		libMesh::ExodusII_IO(*bb_coarsener.get_mesh()).write("prova.e");
	}
}
