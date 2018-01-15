#include "utopia_libmesh.hpp"
#include "utopia_MSHReaderTest.hpp"
#include "utopia_MSHMeshReader.hpp"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"

namespace utopia {
	void test_msh_reader(libMesh::LibMeshInit &init)
	{
		std::cout << "[test_msh_reader]" << std::endl;
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());	

		MeshReaderOpts opts;
		opts.append_mode = true;
		MSHMeshReader reader;

		opts.subdomain_id = 1;
		if(!reader.read("/Users/zulianp/Desktop/algo4u/wearsim/fem.msh", *mesh, opts)) {
			assert(false);
		}

		opts.subdomain_id = 2;
		// if(!reader.read("/Users/zulianp/Desktop/algo4u/wearsim/tibia_insert.msh", *mesh, opts)) {
		if(!reader.read("/Users/zulianp/Desktop/algo4u/wearsim/tibia.msh", *mesh, opts)) {
			assert(false);
		}

		mesh->prepare_for_use();

		libMesh::ExodusII_IO io(*mesh);
		io.write("/Users/zulianp/Desktop/algo4u/wearsim/exodus/wear_geoms.e");
	}
}
