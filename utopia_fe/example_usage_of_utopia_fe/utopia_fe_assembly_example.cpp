
#include "utopia_fe.hpp"
#include <iostream>
#include <memory>

int main(int argc, char *argv[]) {
	using namespace libMesh;
	using namespace std;
	using namespace utopia;

	LibMeshInit init(argc, argv);

	const int n_master = 10;
	const int n_slave = 20;

	auto mesh_master = make_shared<Mesh>(init.comm());
	MeshTools::Generation::build_square (*mesh_master,
		n_master, n_master,
		0, 1,
		0, 1,
		QUAD8);

	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	auto mesh_slave = make_shared<Mesh>(init.comm());
	MeshTools::Generation::build_square (*mesh_slave,
		n_slave, n_slave,
		0.3, 0.8,
		0.3, 0.8,
		QUAD8);

	LibMeshFEContext<LinearImplicitSystem> context_master(mesh_master);
	auto space_master = fe_space(LAGRANGE, Order(1), context_master);

	context_master.equation_systems.init();

	LibMeshFEContext<LinearImplicitSystem> context_slave(mesh_slave);
	auto space_slave = fe_space(LAGRANGE, Order(1), context_slave);

	context_slave.equation_systems.init();

	DSMatrixd B;

	moonolith::Communicator comm(init.comm().get());
	assemble_volume_transfer(comm,
		mesh_master,
		mesh_slave,
		make_ref(space_master.dof_map()),
		make_ref(space_slave.dof_map()),
		0,
		0,
		true,
		1, 
		B);



	return EXIT_SUCCESS;
}
