#ifndef UTOPIA_UI_MORPH_HPP
#define UTOPIA_UI_MORPH_HPP

#include "moonolith_vector.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_ui.hpp"

#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"

#include <set>

namespace utopia {

    template <class Mesh>
    class UIMorph {};  // final : public Configurable { };

    template <>
    class UIMorph<libMesh::DistributedMesh> : public Configurable {
    public:
        UIMorph() : type(""), radius(1.), block(-1) {}

        bool is_valid() const { return !type.empty(); }

        bool apply(libMesh::DistributedMesh &mesh) const {
            if (!is_valid()) return false;

            // std::cout << "morphing " << std::endl;
      #if LIBMESH_VERSION_LESS_THAN(1, 4, 0)
            // Old version
          return false;
      #else
            if (block == -1) {
                auto boundary_node_ids = libMesh::MeshTools::find_boundary_nodes(mesh);

                const libMesh::Point p;
                for (unsigned int n = 0; n < mesh.max_node_id(); n++) {
                    if (boundary_node_ids.count(n)) {
                        auto &node = mesh.node_ref(n);
                        // node *= radius/node.norm();

                        morph_point(node);

                        // std::cout << node(0) << ", " << node(1) << std::endl;
                    }
                }

            } else {
                std::set<libMesh::dof_id_type> morphed;

                for (auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
                    auto &e = **it;

                    if (e.subdomain_id() == block) {
                        if (!e.on_boundary()) continue;

                        int n_sides = e.n_sides();
                        for (int side_num = 0; side_num < n_sides; ++side_num) {
                            if (e.neighbor_ptr(side_num) == nullptr) {
                                auto side_ptr = e.build_side_ptr(side_num);

                                int n_nodes = side_ptr->n_nodes();
                                for (int i = 0; i < n_nodes; ++i) {
                                    auto &node = side_ptr->node_ref(i);
                                    auto id = node.id();

                                    auto m_it = morphed.find(id);

                                    if (m_it == morphed.end()) {
                                        // node *= radius/node.norm();
                                        morph_point(node);
                                        morphed.insert(id);
                                    }
                                }
                            }
                        }
                    }
                }
            }
           #endif
            return true;
        }

        inline void morph_point(libMesh::Point &p) const {
            moonolith::Vector<double, 3> v;
            v.x = p(0);
            v.y = p(1);
            v.z = p(2);

            v = v - center;
            v *= radius / length(v);

            v = center + v;

            p(0) = v.x;
            p(1) = v.y;
            p(2) = v.z;
        }

        void read(Input &is) override {
            // FIXME
            is.get("type", type);
            is.get("radius", radius);
            is.get("block", block);

            is.get("x", center.x);
            is.get("y", center.y);
            is.get("z", center.z);
        }

    private:
        std::string type;
        double radius;
        int block;
        moonolith::Vector<double, 3> center;
    };

}  // namespace utopia

#endif  // UTOPIA_UI_MORPH_HPP
