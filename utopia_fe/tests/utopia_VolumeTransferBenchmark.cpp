#include "utopia_VolumeTransferBenchmark.hpp"
#include "utopia.hpp"

#include "utopia_fe_core.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_assemble_volume_transfer_r.hpp"
#include "utopia_LibMeshBackend.hpp"

#include "moonolith_communicator.hpp"
#include "moonolith_synched_describable.hpp"

#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_refinement.h"

#include <memory>
#include <sstream>

using namespace libMesh;
using std::shared_ptr;
using std::make_shared;

namespace utopia {


	static void refine_mesh(MeshBase &mesh, const int n_refine)
	{
		MeshRefinement mesh_refinement(mesh);
		mesh_refinement.make_flags_parallel_consistent();
		mesh_refinement.uniformly_refine(n_refine);
	}

	void run_experiment(
		LibMeshInit &init, 
		const int n_master,
		const int n_slave,
		const int n_refine_master = 0,
		const int n_refine_slave = 0)
	{
		typedef utopia::LibMeshFEContext<libMesh::LinearImplicitSystem> FEContextT;

		moonolith::Communicator comm(init.comm().get());
		moonolith::root_describe("n procs: " + std::to_string(comm.size()) + "\ncreating fe spaces...", comm, std::cout);
		Chrono c;
		c.start();

		Order order_elem_fine     = FIRST;
		Order order_elem_coarse   = FIRST;
		Order order_of_quadrature = Order(int(order_elem_fine) + int(order_elem_coarse));

		auto slave_mesh = make_shared<DistributedMesh>(init.comm());
		MeshTools::Generation::build_cube(*slave_mesh,
			n_slave, n_slave, n_slave,
			-0.9, 0.9,
			-0.9, 0.9,
			-0.9, 0.9,
			TET4);
			// HEX8);

		auto master_mesh = make_shared<DistributedMesh>(init.comm());
		MeshTools::Generation::build_cube(*master_mesh,
			n_master, n_master, n_master,
			-1., 1.,
			-1., 1.,
			-1., 1.,
			TET4);
			// HEX8);

		if(n_refine_slave) {
			refine_mesh(*slave_mesh, n_refine_slave);
		}

		if(n_refine_master) {
			refine_mesh(*master_mesh, n_refine_master);
		}

		auto slave_context  = make_shared<FEContextT>(slave_mesh);
		auto master_context = make_shared<FEContextT>(master_mesh);

		auto slave_space  = fe_space(LAGRANGE, order_elem_fine, *slave_context);
		auto master_space = fe_space(LAGRANGE, order_elem_coarse, *master_context);

		slave_context->equation_systems.init();
		master_context->equation_systems.init();

		comm.barrier();
		c.stop();
		moonolith::root_describe(c, comm, std::cout);

		std::stringstream ss;
		ss << "experiment: [n_elem: master " << master_mesh->n_elem() << ", slave " << slave_mesh->n_elem() << "]";
		moonolith::root_describe(ss.str(), comm, std::cout);

		c.start();

		DSMatrixd B;
		if(!assemble_volume_transfer(
			comm,
			master_mesh,
			slave_mesh,
			utopia::make_ref(master_space.dof_map()),
			utopia::make_ref(slave_space.dof_map()),
			0,
			0,
			true,
			1,
			B)) 
		{
			assert(false && "Should never get here!");
		}

		const double vol = sum(B);
		const double expected_vol = 1.8*1.8*1.8;
		std::cout << "diff vol: " << std::abs(vol - expected_vol) << std::endl;
		assert(utopia::approxeq(vol, expected_vol, 1e-8));

		comm.barrier();
		c.stop();

		moonolith::root_describe(c, comm, std::cout);
	}

	void run_volume_transfer_benchmark(LibMeshInit &init) 
	{	
		std::vector<std::pair<int,int> > resolutions = 
		{
			{5,  6},
			{10, 10},
			{8,  20},
			{30, 24},
			{10, 44}
		};

		for(auto r : resolutions) {
			run_experiment(init, r.first, r.second, true);
		}
	}

	void run_weak_scaling_benchmark(LibMeshInit &init)
	{
		using std::max;
		const int n_master = max(1., round(pow(mpi_world_size() * 216, 1./3)));
		const int n_slave  = max(1., round(pow(mpi_world_size() * 343, 1./3)));
		run_experiment(init, n_master, n_slave, 1, 1);
	}
}
