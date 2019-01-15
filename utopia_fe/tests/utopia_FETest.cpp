#include "utopia_FETest.hpp"
#include "utopia_BoundaryMeshTest.hpp"

namespace utopia {
	
	void run_all_tests(libMesh::LibMeshInit &init)
	{
		run_boundary_mesh_test(init);
	}

	
	void run_test(const std::string &test_name, libMesh::LibMeshInit &init)
	{
		
	}

}