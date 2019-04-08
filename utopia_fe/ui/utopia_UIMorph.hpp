#ifndef UTOPIA_UI_MORPH_HPP
#define UTOPIA_UI_MORPH_HPP

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_modification.h"

namespace utopia {

    template<class Mesh>
    class UIMorph {};// final : public Configurable { };

    template<>
    class UIMorph<libMesh::DistributedMesh> : public Configurable {
    public:

        UIMorph() : type(""), radius(1.) {}

        bool is_valid() const { return !type.empty(); }

        bool apply(libMesh::DistributedMesh &mesh) const
        {
            if(!is_valid()) return false;

            // std::cout << "morphing " << std::endl;

            auto boundary_node_ids = libMesh::MeshTools::find_boundary_nodes(mesh);

            const libMesh::Point p;
            for (unsigned int n = 0; n < mesh.max_node_id(); n++) {
                if(boundary_node_ids.count(n)) {
                    auto &node = mesh.node_ref(n);
                    node *= radius/node.norm();

                    // std::cout << node(0) << ", " << node(1) << std::endl;
                }
            }

            return true;
        }

        void read(Input &is) override {
                //FIXME
            is.get("type", type);
            is.get("radius", radius);
        }

    private:
        std::string type;
        double radius;
    };

}

#endif //UTOPIA_UI_MORPH_HPP
