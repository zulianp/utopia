#include "utopia_DescribeMeshApp.hpp"

#include "utopia_LibMeshBackend.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_refinement.h"
#include "utopia.hpp"
#include "utopia_libmesh_DiegoMesh.hpp"
#include "utopia_ui.hpp"

namespace utopia {

    DescribeMeshApp::DescribeMeshApp() {}

    DescribeMeshApp::~DescribeMeshApp() {}

    void DescribeMeshApp::run(Input &in) {
        UIMesh<libMesh::DistributedMesh> mesh_reader(comm());
        in.get("mesh", mesh_reader);

        auto &mesh = mesh_reader.mesh();

        std::size_t elem_idx = 0;
        for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
            const std::size_t n_sides = elem_ptr->n_sides();

            for (std::size_t i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                ///////////////////////////// MESH ////////////////////////////////
                auto side_ptr = elem_ptr->build_side_ptr(i);

                int id = utopia::boundary_id(mesh.get_boundary_info(), elem_ptr, i);

                std::cout << id << std::endl;
            }
        }

        // write("desc.e", mesh);

        libMesh::ExodusII_IO(mesh).write("desc.e");
    }
}  // namespace utopia
