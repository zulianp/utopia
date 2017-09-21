#include "utopia_VolumeTransferBenchmark.hpp"
#include "utopia.hpp"

#include "utopia_fe_core.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_assemble_volume_transfer_r.hpp"
#include "utopia_LibMeshBackend.hpp"

#include "moonolith_communicator.hpp"

#include <libmesh/parallel_mesh.h>
#include <memory>

using namespace libMesh;
using std::shared_ptr;
using std::make_shared;

namespace utopia {
	void run_volume_transfer_benchmark(LibMeshInit &init) {
		typedef utopia::LibMeshFEContext<libMesh::LinearImplicitSystem> FEContextT;

		Order order_elem_fine     = FIRST;
		Order order_elem_coarse   = FIRST;
		Order order_of_quadrature = Order(int(order_elem_fine) + int(order_elem_coarse));

		const std::string data_path = Utopia::Instance().get("data_path");
		const std::string fine_mesh_path   = data_path + "/fine_mesh.e";
		const std::string coarse_mesh_path = data_path + "/coarse_mesh.e";

		auto fine_mesh = make_shared<DistributedMesh>(init.comm());
		fine_mesh->read(fine_mesh_path);

		auto coarse_mesh = make_shared<DistributedMesh>(init.comm());
		coarse_mesh->read(coarse_mesh_path);

		auto fine_context   = make_shared<FEContextT>(fine_mesh);
		auto coarse_context = make_shared<FEContextT>(coarse_mesh);

		auto fine_space   = fe_space(LAGRANGE, order_elem_fine, *fine_context);
		auto coarse_space = fe_space(LAGRANGE, order_elem_coarse, *coarse_context);

		auto u = fe_function(fine_space);
		strong_enforce( boundary_conditions(u == coeff(0.), {1}) );

		fine_context->equation_systems.init();
		coarse_context->equation_systems.init();

		moonolith::Communicator comm(init.comm().get());
		DSMatrixd B;

		if(!assemble_volume_transfer(
			comm,
			coarse_mesh,
			fine_mesh,
			utopia::make_ref(coarse_space.dof_map()),
			utopia::make_ref(fine_space.dof_map()),
			0,
			0,
			true,
			1,
			B)) {

			assert(false);
		}
	}
}
