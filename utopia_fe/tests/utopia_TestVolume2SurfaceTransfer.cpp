#include "utopia_TestVolume2SurfaceTransfer.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"
#include "moonolith_communicator.hpp"
#include "utopia_assemble_volume_transfer.hpp"

typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

namespace utopia {
	void run_volume_to_surface_transfer_test(libMesh::LibMeshInit &init)
	{
		auto n = 1;
		auto elem_type  = libMesh::HEX8;
		auto elem_order = libMesh::FIRST;

		auto vol_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		libMesh::MeshTools::Generation::build_cube(
			*vol_mesh,
			n, n, n,
			0, 1,
			0, 1.,
			0, 1,
			elem_type
		);

		auto surf_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		surf_mesh->read("../data/test/square_with_2_tri.e");

		//equations system
		auto vol_equation_systems = std::make_shared<libMesh::EquationSystems>(*vol_mesh);
		auto &vol_sys = vol_equation_systems->add_system<libMesh::LinearImplicitSystem>("vol_sys");

		auto surf_equation_systems = std::make_shared<libMesh::EquationSystems>(*surf_mesh);
		auto &surf_sys = surf_equation_systems->add_system<libMesh::LinearImplicitSystem>("surf_sys");

		//scalar function space
		auto V_vol  = FunctionSpaceT(vol_equation_systems, libMesh::LAGRANGE, elem_order,   "u_vol");
		auto V_surf = FunctionSpaceT(surf_equation_systems, libMesh::LAGRANGE, elem_order, "u_surf");

		V_vol.initialize();
		V_surf.initialize();

		DSMatrixd B;
		moonolith::Communicator comm(init.comm().get());
		if(assemble_volume_transfer(
		    comm,
		    vol_mesh,
		    surf_mesh,
		    make_ref(V_vol.dof_map()),
		    make_ref(V_surf.dof_map()),
		    0,
		    0,
		    false, 
		    1,
		    B))
		{
			write("B.m", B);

			DVectord v_vol = local_values(V_vol.dof_map().n_local_dofs(), 1.);
			DVectord v_surf = B * v_vol;

			disp(v_vol);
			disp(v_surf);

			// Nemesis_IO vol_IO(*vol_mesh);
			// vol_IO.write_equation_systems ("surf2vol_vol.e", vol_sys);

			// Nemesis_IO surf_IO(*surf_mesh);
			// surf_IO.write_equation_systems ("surf2vol_surf.e", surf_sys);

		} else {
			assert(false);
		}
	}
}