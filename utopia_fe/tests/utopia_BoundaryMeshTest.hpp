#ifndef UTOPIA_BOUNDARY_MESH_TEST_HPP
#define UTOPIA_BOUNDARY_MESH_TEST_HPP

#include "utopia_fe_base.hpp"



namespace libMesh {
	class LibMeshInit;
}


namespace utopia {
	void run_boundary_mesh_test(libMesh::LibMeshInit &init);
}



#endif //UTOPIA_BOUNDARY_MESH_TEST_HPP
