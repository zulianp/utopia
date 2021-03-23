// #include "utopia_DescribeMeshApp.hpp"
#include "utopia_Main.hpp"

#include "utopia_LibMeshBackend.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_refinement.h"
#include "utopia.hpp"
#include "utopia_libmesh_DiegoMesh.hpp"
#include "utopia_ui.hpp"

#include "utopia_libmesh_Mesh.hpp"

namespace utopia {

    void describe_mesh(Input &in) {
        utopia::libmesh::Mesh mesh_wrapper;
        in.get("mesh", mesh_wrapper);

        if (mesh_wrapper.empty()) {
            utopia::err() << "[Error] something went wrong reading the mesh\n";
            return;
        }

        auto &mesh = mesh_wrapper.raw_type();

        std::size_t elem_idx = 0;
        for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
            const std::size_t n_sides = elem_ptr->n_sides();

            for (std::size_t i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                ///////////////////////////// MESH ////////////////////////////////
                auto side_ptr = elem_ptr->build_side_ptr(i);

                auto &bi = mesh.get_boundary_info();
                int id = utopia::boundary_id(bi, elem_ptr, i);

                auto name = bi.get_sideset_name(id);

                std::cout << id << " " << name << " " << bi.get_id_by_name(name) << std::endl;
            }
        }

        // write("desc.e", mesh);

        libMesh::ExodusII_IO(mesh).write("desc.e");
    }

    UTOPIA_REGISTER_APP(describe_mesh);
}  // namespace utopia
