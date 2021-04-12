#include "utopia_MSHReaderTest.hpp"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "utopia_MSHMeshReader.hpp"
#include "utopia_UGMeshReader.hpp"
#include "utopia_libmesh_old.hpp"

namespace utopia {
    void MSHReaderTest::run(Input &in) {
        std::cout << "[test_msh_reader]" << std::endl;
        auto mesh = std::make_shared<libMesh::DistributedMesh>(this->comm());

        // MeshReaderOpts opts;
        // opts.append_mode = true;
        // MSHMeshReader reader;

        // opts.subdomain_id = 1;
        // if(!reader.read("/Users/zulianp/Desktop/algo4u/wearsim/fem.msh", *mesh, opts)) {
        // 	assert(false);
        // }

        // opts.subdomain_id = 2;
        // // if(!reader.read("/Users/zulianp/Desktop/algo4u/wearsim/tibia_insert.msh", *mesh, opts)) {
        // if(!reader.read("/Users/zulianp/Desktop/algo4u/wearsim/tibia.msh", *mesh, opts)) {
        // 	assert(false);
        //  return;
        // }

        // mesh->prepare_for_use();

        UGXMeshReader reader;
        if (!reader.read("/Users/zulianp/Desktop/algo4u/wearsim/promesh/tibia_insert_wp.ugx", *mesh)) {
            utopia_test_assert(false);
            return;
        }

        libMesh::ExodusII_IO io(*mesh);
        io.write("/Users/zulianp/Desktop/algo4u/wearsim/exodus/tibia_insert.e");
    }
}  // namespace utopia
