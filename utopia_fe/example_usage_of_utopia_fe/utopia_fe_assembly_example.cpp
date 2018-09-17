
#include "utopia_fe.hpp"
#include <iostream>
#include <memory>

#include "utopia_libmesh.hpp"

int main(int argc, char *argv[]) {
	using namespace libMesh;
	using namespace std;
	using namespace utopia;

	LibMeshInit init(argc, argv);

	const int n_master = 10;
	const int n_slave = 20;

	auto mesh_master = make_shared<DistributedMesh>(init.comm());
	MeshTools::Generation::build_square (*mesh_master,
		n_master, n_master,
		0, 1,
		0, 1,
		QUAD8);

	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	auto mesh_slave = make_shared<DistributedMesh>(init.comm());
	MeshTools::Generation::build_square (*mesh_slave,
		n_slave, n_slave,
		0.3, 0.8,
		0.3, 0.8,
		QUAD8);

	UVector d = local_zeros(10);

	disp(d);


	return EXIT_SUCCESS;
}
