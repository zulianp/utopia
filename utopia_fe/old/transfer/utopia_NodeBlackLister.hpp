#ifndef UTOPIA_NODE_BLACK_LISTER_HPP
#define UTOPIA_NODE_BLACK_LISTER_HPP

#include "utopia_libmesh.hpp"
#include "utopia_ui.hpp"

#include <set>
#include <vector>

#include <libmesh/elem.h>
#include <libmesh/mesh_base.h>

namespace utopia {

    class NodeSetBlackLister final : public Configurable {
    public:
        void read(Input &in) override;
        bool is_black_listed(const int node_set_id) const;
        bool is_black_listed(const libMesh::MeshBase &mesh, const libMesh::Elem &e) const;
        bool has_black_listed_nodes(const libMesh::MeshBase &mesh) const;

    public:
        std::set<int> node_sets_;
    };

    class ElementBlackList : public Configurable {
    public:
        void read(Input &in) override;
        void init(const libMesh::MeshBase &mesh);
        bool is_black_listed(const SizeType element_id);
        bool is_black_listed(const libMesh::Elem &e);
        ElementBlackList(const bool is_boundary_mesh);

    private:
        NodeSetBlackLister black_lister_;
        bool is_boundary_mesh_;
        std::set<SizeType> black_listed_;

        void init(const NodeSetBlackLister &black_lister, const libMesh::MeshBase &mesh);
    };

}  // namespace utopia

#endif  // UTOPIA_NODE_BLACK_LISTER_HPP
